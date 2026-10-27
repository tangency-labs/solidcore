//! Generates `UpstreamFixes` (`src/kernel.rs`) from one file per fix under `fixes/`.
//!
//! Each opt-in fix of an upstream defect is `fixes/<field_name>.toml`, so pull requests
//! that add different fixes add different files instead of all editing one struct. The
//! files use a small subset of TOML, exactly these three keys in any order:
//!
//! ```toml
//! records = ["WR-300"]        # the plan/upstream-issues records, may be empty
//! issues = [1023]             # the tangency-labs/solidcore issues, may be empty
//! doc = '''
//! The field's doc comment, one or more lines, written as rustdoc (WR-300, #1023).
//! '''
//! ```
//!
//! At least one record or issue is required, and each must appear in `doc` (as `WR-300`
//! and `#1023`). Comment lines (`#`) and blank lines are allowed between keys. The file
//! stem is the field name: lower-case snake case and not a Rust keyword.
//!
//! The output, `$OUT_DIR/upstream_fixes.rs`, defines the macro `define_upstream_fixes!`,
//! which `src/kernel.rs` invokes with the struct's doc comment, attributes and name. The
//! files are read in sorted order, so the output is deterministic.

use std::fmt::Write as _;
use std::path::Path;

/// One parsed `fixes/<name>.toml`.
struct Fix {
    name: String,
    records: Vec<String>,
    issues: Vec<u32>,
    doc: Vec<String>,
}

const KEYWORDS: &[&str] = &[
    "as", "async", "await", "break", "const", "continue", "crate", "dyn", "else", "enum", "extern",
    "false", "fn", "for", "gen", "if", "impl", "in", "let", "loop", "match", "mod", "move", "mut",
    "pub", "ref", "return", "self", "static", "struct", "super", "trait", "true", "try", "type",
    "unsafe", "use", "where", "while", "abstract", "become", "box", "do", "final", "macro",
    "override", "priv", "typeof", "unsized", "virtual", "yield",
];

fn main() -> Result<(), String> {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=fixes");
    let dir = Path::new("fixes");
    let mut paths = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))? {
        let path = entry.map_err(|e| format!("{}: {e}", dir.display()))?.path();
        paths.push(path);
    }
    paths.sort();
    let mut fixes: Vec<Fix> = Vec::new();
    for path in &paths {
        let shown = path.display().to_string();
        let name = match (path.file_stem(), path.extension()) {
            (Some(stem), Some(ext)) if ext == "toml" => stem
                .to_str()
                .ok_or_else(|| format!("{shown}: file name is not UTF-8"))?
                .to_owned(),
            _ => {
                return Err(format!(
                    "{shown}: only `<field_name>.toml` files belong in fixes/"
                ))
            }
        };
        check_ident(&name).map_err(|e| format!("{shown}: {e}"))?;
        if fixes.iter().any(|f| f.name == name) {
            return Err(format!("{shown}: duplicate fix `{name}`"));
        }
        let text = std::fs::read_to_string(path).map_err(|e| format!("{shown}: {e}"))?;
        fixes.push(parse(name, &text).map_err(|e| format!("{shown}: {e}"))?);
    }
    if fixes.is_empty() {
        return Err("fixes/: no fixes found".to_owned());
    }
    let out = std::env::var("OUT_DIR").map_err(|e| format!("OUT_DIR: {e}"))?;
    let out = Path::new(&out).join("upstream_fixes.rs");
    std::fs::write(&out, generate(&fixes)).map_err(|e| format!("{}: {e}", out.display()))
}

/// Checks that `name` is a lower-case snake-case identifier and not a keyword.
fn check_ident(name: &str) -> Result<(), String> {
    let mut chars = name.chars();
    let first_ok = chars
        .next()
        .is_some_and(|c| c.is_ascii_lowercase() || c == '_');
    let rest_ok = chars.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_');
    if !first_ok || !rest_ok || name == "_" {
        return Err(format!(
            "`{name}` is not a lower-case snake-case field name"
        ));
    }
    if KEYWORDS.contains(&name) {
        return Err(format!("`{name}` is a Rust keyword"));
    }
    Ok(())
}

/// Parses one fix file (the TOML subset in the module docs).
fn parse(name: String, text: &str) -> Result<Fix, String> {
    let mut records = None;
    let mut issues = None;
    let mut doc = None;
    let mut lines = text.lines().enumerate();
    while let Some((n, line)) = lines.next() {
        let at = n + 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let (key, value) = trimmed
            .split_once('=')
            .ok_or_else(|| format!("line {at}: expected `key = value`"))?;
        let (key, value) = (key.trim(), value.trim());
        match key {
            "records" => {
                let mut list = Vec::new();
                for item in array(value).map_err(|e| format!("line {at}: {e}"))? {
                    let record = item
                        .strip_prefix('"')
                        .and_then(|s| s.strip_suffix('"'))
                        .filter(|s| is_record(s))
                        .ok_or_else(|| format!("line {at}: record `{item}` is not \"WR-<n>\""))?;
                    list.push(record.to_owned());
                }
                set_once(&mut records, list, key, at)?;
            }
            "issues" => {
                let mut list = Vec::new();
                for item in array(value).map_err(|e| format!("line {at}: {e}"))? {
                    let issue = item
                        .parse::<u32>()
                        .ok()
                        .filter(|&i| i > 0 && !item.starts_with('0'))
                        .ok_or_else(|| format!("line {at}: issue `{item}` is not a number"))?;
                    list.push(issue);
                }
                set_once(&mut issues, list, key, at)?;
            }
            "doc" => {
                if value != "'''" {
                    return Err(format!(
                        "line {at}: `doc` must open with `'''` and continue on the next line"
                    ));
                }
                let mut body = Vec::new();
                loop {
                    let (_, l) = lines
                        .next()
                        .ok_or_else(|| format!("line {at}: `doc` has no closing `'''`"))?;
                    if l.trim_end() == "'''" {
                        break;
                    }
                    if l.contains("'''") {
                        return Err(format!(
                            "line {at}: `doc`'s closing `'''` must be on its own line"
                        ));
                    }
                    body.push(l.trim_end().to_owned());
                }
                if body.iter().all(|l| l.trim().is_empty()) {
                    return Err(format!("line {at}: `doc` is empty"));
                }
                set_once(&mut doc, body, key, at)?;
            }
            _ => return Err(format!("line {at}: unknown key `{key}`")),
        }
    }
    let records = records.ok_or("missing `records`")?;
    let issues = issues.ok_or("missing `issues`")?;
    let doc = doc.ok_or("missing `doc`")?;
    if records.is_empty() && issues.is_empty() {
        return Err("name at least one record or issue".to_owned());
    }
    let joined = doc.join("\n");
    for r in &records {
        if !mentions(&joined, r) {
            return Err(format!("`doc` does not mention record {r}"));
        }
    }
    for i in &issues {
        if !mentions(&joined, &format!("#{i}")) {
            return Err(format!("`doc` does not mention issue #{i}"));
        }
    }
    Ok(Fix {
        name,
        records,
        issues,
        doc,
    })
}

fn set_once<T>(slot: &mut Option<T>, value: T, key: &str, at: usize) -> Result<(), String> {
    if slot.replace(value).is_some() {
        return Err(format!("line {at}: duplicate key `{key}`"));
    }
    Ok(())
}

/// The trimmed items of a one-line `[a, b]` array.
fn array(value: &str) -> Result<Vec<&str>, String> {
    let inner = value
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .ok_or("expected a one-line `[...]` array")?
        .trim();
    if inner.is_empty() {
        return Ok(Vec::new());
    }
    inner
        .split(',')
        .map(|s| {
            let s = s.trim();
            if s.is_empty() {
                Err("empty array item".to_owned())
            } else {
                Ok(s)
            }
        })
        .collect()
}

fn is_record(s: &str) -> bool {
    s.strip_prefix("WR-")
        .is_some_and(|n| !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()))
}

/// Whether `text` contains `token` not followed by another digit.
fn mentions(text: &str, token: &str) -> bool {
    text.match_indices(token).any(|(i, _)| {
        !text
            .get(i + token.len()..)
            .and_then(|rest| rest.chars().next())
            .is_some_and(|c| c.is_ascii_digit())
    })
}

/// The generated `define_upstream_fixes!` macro.
fn generate(fixes: &[Fix]) -> String {
    let mut s = String::new();
    s.push_str("// @generated by crates/solidcore/build.rs from crates/solidcore/fixes/*.toml.\n");
    s.push_str("// Add a fix by adding a file there; do not edit this output.\n\n");
    s.push_str("macro_rules! define_upstream_fixes {\n");
    s.push_str("    ($(#[$meta:meta])* $vis:vis struct $name:ident;) => {\n");
    s.push_str("        $(#[$meta])*\n        $vis struct $name {\n");
    for f in fixes {
        for line in &f.doc {
            let line = if line.is_empty() {
                String::new()
            } else {
                format!(" {line}")
            };
            let _ = writeln!(s, "            #[doc = {line:?}]");
        }
        let _ = writeln!(s, "            pub {}: bool,", f.name);
    }
    s.push_str("        }\n\n        impl $name {\n");
    s.push_str(
        "            /// Every opt-in fix, in field-name order: its field name, records and\n",
    );
    s.push_str("            /// issues (one file each under `crates/solidcore/fixes/`).\n");
    s.push_str("            pub const ALL: &'static [FixInfo] = &[\n");
    for f in fixes {
        let records: Vec<String> = f.records.iter().map(|r| format!("{r:?}")).collect();
        let issues: Vec<String> = f.issues.iter().map(u32::to_string).collect();
        let _ = writeln!(
            s,
            "                FixInfo {{ name: {:?}, records: &[{}], issues: &[{}] }},",
            f.name,
            records.join(", "),
            issues.join(", ")
        );
    }
    s.push_str("            ];\n\n");
    s.push_str("            /// Every fix on.\n");
    s.push_str("            #[must_use]\n");
    s.push_str("            pub fn all_on() -> Self {\n                Self {\n");
    for f in fixes {
        let _ = writeln!(s, "                    {}: true,", f.name);
    }
    s.push_str("                }\n            }\n\n");
    s.push_str("            /// The field called `name`, or `None` for an unknown name.\n");
    s.push_str("            pub fn get(&self, name: &str) -> Option<bool> {\n");
    s.push_str("                match name {\n");
    for f in fixes {
        let _ = writeln!(
            s,
            "                    {:?} => Some(self.{}),",
            f.name, f.name
        );
    }
    s.push_str("                    _ => None,\n                }\n            }\n\n");
    s.push_str(
        "            /// The field called `name`, to set it, or `None` for an unknown name.\n",
    );
    s.push_str("            pub fn get_mut(&mut self, name: &str) -> Option<&mut bool> {\n");
    s.push_str("                match name {\n");
    for f in fixes {
        let _ = writeln!(
            s,
            "                    {:?} => Some(&mut self.{}),",
            f.name, f.name
        );
    }
    s.push_str(
        "                    _ => None,\n                }\n            }\n        }\n    };\n}\n",
    );
    s
}
