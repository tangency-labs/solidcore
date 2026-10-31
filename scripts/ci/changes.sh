#!/usr/bin/env bash
# Prints the CI `changes` job outputs (code=, full=) for .github/workflows/ci.yml.
# code=false only for a pull_request or merge_group whose diff against its base touches
# nothing outside DOCS_ONLY (prose read by no build, test or difftest) and no upstream-issue
# record a known mismatch's `cpp_defect` names. Anything else, and
# any doubt (no base, fetch failure), is code=true. full=true when code=true and the event
# is not a push to main (pushes run only the jobs that publish main's artifacts).
# Usage: changes.sh <event> [base sha]
set -euo pipefail
event=${1:?event}; base=${2:-}
# Paths no cargo build/test, difftest run or packaging reads. Not plan/units/ (difftest reads
# unit statuses), not crates/** READMEs (cargo package reads them), not scripts/ or harness/.
# plan/upstream-issues/ records a known mismatch names in `cpp_defect` are excepted below.
DOCS_ONLY='^(plan/findings/.*|plan/upstream-issues/.*|plan/[^/]*\.md|site/.*|docs/.*|assets/.*|[^/]*\.md)$'

# The ids known mismatches name in `cpp_defect`. difftest reads each one's record,
# plan/upstream-issues/<id>--<slug>.md, but not its addenda (crates/difftest/src/known.rs,
# check_cpp_defect).
cpp_defect_ids() {
  sed -n 's/^[[:space:]]*cpp_defect[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' \
    crates/difftest/known_mismatches/*.toml 2>/dev/null | sort -u || true
}
code=true
if [ "$event" = pull_request ] || [ "$event" = merge_group ]; then
  if [ -n "$base" ] && { git cat-file -e "$base^{commit}" 2>/dev/null || git fetch -q --depth=1 origin "$base"; }; then
    # --no-renames: a rename lists both its old and new path, so moving code into a docs path
    # still counts as a code change. Quoted (unusual) paths never match, which fails closed.
    files=$(git diff --no-renames --name-only "$base" HEAD)
    if [ -n "$files" ] && ! grep -qvE "$DOCS_ONLY" <<<"$files"; then
      code=false
    fi
    # Editing, renaming or deleting such a record, or adding a second one with its id, can
    # break the difftest gate, so it is a code change.
    if [ "$code" = false ]; then
      for id in $(cpp_defect_ids); do
        while IFS= read -r f; do
          case "$f" in
            *--addendum-*) ;;
            "plan/upstream-issues/$id--"*) code=true ;;
          esac
        done <<<"$files"
      done
    fi
    echo "changed files: $(wc -l <<<"$files"), code=$code" >&2
  fi
fi
full=$code
[ "$event" = push ] && full=false
echo "code=$code"
echo "full=$full"
