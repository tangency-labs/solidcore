//! Where solidcore's golden data lives: the goldens (`expected/<unit>.jsonl`) and the per-unit
//! fixtures (`fixtures/<unit>/`), kept in tangency-labs/solidcore-goldens and checked out at the
//! sha in the repository's `goldens.pin`.
//!
//! The data root is `$SOLIDCORE_GOLDENS_DIR`, else `.goldens/` at the repository root, which
//! `cargo xtask goldens-fetch` fills. Every reader goes through [`root`]: a root that is missing,
//! lacks `expected/` or `fixtures/`, or is a checkout at another sha than the pin is an error
//! naming the fetch command. Nothing here skips: a test that needs data and has none fails.
//!
//! `fixtures/solidcraft/` (the SolidCraft case registry) is not golden data and stays in the
//! repository.

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

/// The environment variable naming the data root.
pub const ENV: &str = "SOLIDCORE_GOLDENS_DIR";

/// The data root's default name, at the repository root.
pub const DEFAULT_DIR: &str = ".goldens";

/// The file at the repository root holding the pinned sha.
pub const PIN_FILE: &str = "goldens.pin";

/// The command that fetches the data at the pin.
pub const FETCH: &str = "cargo xtask goldens-fetch";

/// Why the golden data cannot be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Missing {
    message: String,
}

impl Missing {
    fn new(message: String) -> Self {
        Self { message }
    }
}

impl fmt::Display for Missing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Missing {}

/// The repository root (two levels above this crate's manifest).
#[must_use]
pub fn repo_root() -> PathBuf {
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    here.parent()
        .and_then(Path::parent)
        .map_or_else(|| here.to_path_buf(), Path::to_path_buf)
}

/// The pinned sha, from `goldens.pin` at the repository root (its first word).
///
/// # Errors
/// [`Missing`] when the file cannot be read or holds no 40-digit hex sha.
pub fn pin() -> Result<String, Missing> {
    let path = repo_root().join(PIN_FILE);
    let text =
        fs::read_to_string(&path).map_err(|e| Missing::new(format!("{}: {e}", path.display())))?;
    let sha = text.split_whitespace().next().unwrap_or_default();
    if sha.len() == 40 && sha.bytes().all(|b| b.is_ascii_hexdigit()) {
        Ok(sha.to_ascii_lowercase())
    } else {
        Err(Missing::new(format!(
            "{}: not a full sha: `{sha}`",
            path.display()
        )))
    }
}

/// The sha a checkout at `dir` stands at (a detached `HEAD`), if it is one.
fn checked_out(dir: &Path) -> Option<String> {
    let head = fs::read_to_string(dir.join(".git").join("HEAD")).ok()?;
    let sha = head.trim();
    (sha.len() == 40 && sha.bytes().all(|b| b.is_ascii_hexdigit()))
        .then(|| sha.to_ascii_lowercase())
}

/// The data root, checked: `$SOLIDCORE_GOLDENS_DIR` or `.goldens/` at the repository root,
/// holding `expected/` and `fixtures/`, and, when it is a git checkout, at the pinned sha.
///
/// # Errors
/// [`Missing`], naming [`FETCH`], for a root that is absent, incomplete or at another sha.
pub fn root() -> Result<PathBuf, Missing> {
    let (dir, from) = match std::env::var_os(ENV) {
        Some(d) if !d.is_empty() => (PathBuf::from(d), format!("${ENV}")),
        _ => (repo_root().join(DEFAULT_DIR), DEFAULT_DIR.to_string()),
    };
    for sub in ["expected", "fixtures"] {
        if !dir.join(sub).is_dir() {
            return Err(Missing::new(format!(
                "golden data not found: {} ({from}) has no `{sub}/`; run `{FETCH}` \
                 (or set ${ENV} to a checkout of tangency-labs/solidcore-goldens)",
                dir.display()
            )));
        }
    }
    if let Some(at) = checked_out(&dir) {
        let want = pin()?;
        if at != want {
            return Err(Missing::new(format!(
                "golden data at {} is {at}, but {PIN_FILE} pins {want}; run `{FETCH}`",
                dir.display()
            )));
        }
    }
    Ok(dir)
}

/// The goldens directory, `<root>/expected`.
///
/// # Errors
/// Those of [`root`].
pub fn expected_dir() -> Result<PathBuf, Missing> {
    Ok(root()?.join("expected"))
}

/// The per-unit fixtures directory, `<root>/fixtures`.
///
/// # Errors
/// Those of [`root`].
pub fn fixtures_dir() -> Result<PathBuf, Missing> {
    Ok(root()?.join("fixtures"))
}

/// The path of `rel` (`fixtures/...` or `expected/...`, no `..`) under the data root.
///
/// # Errors
/// [`Missing`] for another `rel`; those of [`root`].
pub fn path(rel: &str) -> Result<PathBuf, Missing> {
    let ok = (rel.starts_with("fixtures/") || rel.starts_with("expected/"))
        && !rel.split(['/', '\\']).any(|c| c == "..");
    if !ok {
        return Err(Missing::new(format!(
            "`{rel}`: golden data paths are under `fixtures/` or `expected/`, without `..`"
        )));
    }
    Ok(root()?.join(rel))
}

/// The bytes of `rel` under the data root ([`path`]).
///
/// # Errors
/// Those of [`path`]; [`Missing`] naming the file when it cannot be read.
pub fn read(rel: &str) -> Result<Vec<u8>, Missing> {
    let p = path(rel)?;
    fs::read(&p).map_err(|e| Missing::new(format!("{}: {e}", p.display())))
}

/// The text of `rel` under the data root ([`path`]).
///
/// # Errors
/// Those of [`read`]; [`Missing`] for text that is not UTF-8.
pub fn read_to_string(rel: &str) -> Result<String, Missing> {
    let p = path(rel)?;
    fs::read_to_string(&p).map_err(|e| Missing::new(format!("{}: {e}", p.display())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_pin_is_a_full_sha() {
        assert_eq!(pin().unwrap().len(), 40);
    }

    #[test]
    fn paths_outside_the_data_are_refused() {
        for rel in [
            "Cargo.toml",
            "fixtures/../Cargo.toml",
            "upstream/x",
            "fixtures\\..\\x",
        ] {
            let e = path(rel).unwrap_err().to_string();
            assert!(e.contains("golden data paths"), "{rel}: {e}");
        }
    }
}
