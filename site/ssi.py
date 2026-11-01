"""Render measured SSI results; absent or invalid artifacts never become zero scores."""
from datetime import datetime, timezone
from html import escape as esc
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import ssi_report


def load(root, repository):
    try:
        return ssi_report.validate_bundle(root / "reports/ssi", repository, root)
    except (ssi_report.ReportError, OSError) as error:
        print(f"site: SSI report unavailable: {error}", file=sys.stderr)
        return None, None


def engines(report):
    return [engine for engine in ssi_report.ENGINES if engine in report.get("summary", {})]


def summary_table(report):
    names = engines(report)
    summary = report.get("summary", {})
    verdicts = sorted({v for engine in names for v in summary[engine]})
    head = "".join(f"<th>{esc(ssi_report.LABELS[e])}</th>" for e in names)
    rows = "".join("<tr><td>" + esc(v) + "</td>" + "".join(
        f'<td class="num">{summary[e].get(v, 0)}</td>' for e in names) + "</tr>" for v in verdicts)
    return ('<div class="scroll narrow"><table><thead><tr><th>Verdict</th>' + head
            + "</tr></thead><tbody>" + rows + "</tbody></table></div>")


def measurement(provenance, repository, current_upstream):
    if not provenance:
        return ""
    sha, upstream = provenance["source_sha"], provenance["upstream_commit"]
    date = datetime.fromisoformat(provenance["generated_at"].replace("Z", "+00:00"))
    date = date.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    text = (f'<p class="muted">Measured {esc(date)} on '
            f'<a href="https://github.com/{esc(repository)}/commit/{sha}"><code>{sha[:12]}</code></a>, '
            f'SMLib v{esc(provenance["upstream_version"])} (<code>{upstream[:12]}</code>). '
            f'<a href="{esc(provenance["run_url"])}">Benchmark run</a> · '
            '<a href="reports/ssi-corpus.json">Report JSON</a> · '
            '<a href="reports/ssi-provenance.json">Provenance JSON</a>.</p>')
    if upstream != current_upstream:
        text += '<p class="warn">These results use a different SMLib pin from the current site.</p>'
    return text


def unavailable(repository):
    return ('<p class="muted">No SSI report available. The complete benchmark runs nightly '
            f'and <a href="https://github.com/{esc(repository)}/actions/workflows/ssi-corpus.yml">'
            'can be run manually</a>. Results appear after a validated artifact is available.</p>')


def overview(report, provenance, repository, current_upstream):
    text = '<section><h2>Surface/surface intersection corpus</h2>'
    if report is None:
        return text + unavailable(repository) + '<p><a href="ssi.html">All cases</a>.</p></section>'
    parity = report.get("summary", {}).get("parity", {})
    text += (f'<p>{len(report.get("cases", []))} cases scored against known true intersections. '
             f'C++/Rust parity: {parity.get("compared", 0)} compared, '
             f'{parity.get("identical", 0)} agree within 1e-9, '
             f'{parity.get("mismatched", 0)} mismatched. <a href="ssi.html">All cases</a>.</p>')
    return text + measurement(provenance, repository, current_upstream) + summary_table(report) + "</section>"


def body(report, provenance, repository, current_upstream):
    text = "<h1>SSI corpus</h1>"
    if report is None:
        return text + unavailable(repository)
    names = engines(report)
    cases = sorted(report.get("cases", []), key=lambda case: case["id"])
    corpus = report.get("corpus", {})
    text += measurement(provenance, repository, current_upstream)
    text += (f'<p class="muted">{len(cases)} surface/surface intersection cases from '
             f'<a href="{esc(corpus.get("repo", ""))}">ssi-corpus</a> at '
             f'<code>{esc(corpus.get("commit", "")[:12])}</code>. '
             f'Approximation tolerance {esc(str(report.get("approx_tol", "")))}, accuracy tolerance '
             f'{esc(str(report.get("acc_tol", "")))}. C++ and Rust sample each curve at '
             f'{report.get("samples_per_curve", 0)} parameters. '
             'Parity compares C++ with the Rust port; ssi and solidcore-ssi are scored independently.</p>')
    text += summary_table(report)
    classes = sorted({case.get("class", "") for case in cases})
    class_options = "".join(f'<option value="{esc(c)}">{esc(c)}</option>' for c in classes)
    text += ('<div class="filters" data-table="ssi"><label>Search <input type="search" name="q" '
             'placeholder="case id or surface"></label><label>Class <select name="class">'
             '<option value="">all</option>' + class_options + '</select></label>')
    for engine in names:
        options = "".join(f'<option value="{esc(v)}">{esc(v)} ({n})</option>'
                          for v, n in sorted(report["summary"][engine].items()))
        text += (f'<label>{esc(ssi_report.LABELS[engine])} <select name="{engine}">'
                 '<option value="">all</option>' + options + '</select></label>')
    text += '<span class="count" aria-live="polite"></span></div>'
    headings = "".join(f'<th data-sort="text">{esc(ssi_report.LABELS[e])}</th>' for e in names)
    text += ('<div class="scroll"><table id="ssi" class="sortable"><thead><tr>'
             '<th data-sort="text">Case</th><th data-sort="text">Class</th>'
             '<th data-sort="text">Surfaces</th>' + headings
             + '<th data-sort="text">C++/Rust parity</th><th>Flags</th></tr></thead><tbody>')
    for case in cases:
        attributes = "".join(f' data-{e}="{esc(case[e].get("verdict", ""))}"' for e in names)
        row = (f'<tr data-class="{esc(case.get("class", ""))}"{attributes}>'
               f'<td><code>{esc(case["id"])}</code></td><td>{esc(case.get("class", ""))}</td>'
               f'<td>{esc(" × ".join(case.get("surfaces", [])))}</td>')
        flags, cpp_flags = [], case.get("cpp", {}).get("flags", [])
        for engine in names:
            grade = case[engine]
            verdict = grade.get("verdict", "")
            row += f'<td><span class="v v-{esc(verdict.replace(" ", "-"))}">{esc(verdict)}</span>'
            if grade.get("status_name"):
                row += f'<div class="small">{esc(grade["status_name"])}</div>'
            if grade.get("deferred"):
                deferred = grade["deferred"]
                unit = esc(deferred.get("unit", ""))
                row += (f'<div class="small">waits on <a href="units/{unit}.html">{unit}</a>: '
                        f'<code>{esc(deferred.get("symbol", ""))}</code></div>')
            row += "</td>"
            details = grade.get("flags", [])
            if details and (engine == "cpp" or details != cpp_flags):
                flags.extend(f'<li>{esc(ssi_report.LABELS[engine])}: {esc(flag)}</li>' for flag in details)
            if grade.get("detail"):
                flags.append(f'<li>{esc(ssi_report.LABELS[engine])}: {esc(grade["detail"])}</li>')
        parity = case.get("parity")
        agreement = ("agree" if parity["pass"] else f'{parity["mismatches"]} mismatches') if parity else ""
        text += row + f'<td>{esc(agreement)}</td><td class="flags">'
        text += "<ul>" + "".join(flags) + "</ul>" if flags else ""
        text += "</td></tr>"
    return text + "</tbody></table></div>"
