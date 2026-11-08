//! The golden data's commands. `cargo xtask goldens-fetch [--check]`: the golden data (tangency-labs/solidcore-goldens) at
//! the sha in `goldens.pin`, shallow-fetched into the data root ([`goldens_data::root`]'s:
//! `$SOLIDCORE_GOLDENS_DIR`, else `.goldens/`). Idempotent: a checkout already at the pin with
//! no changes is left alone; one at another sha is moved to the pin; one with changes (a bless
//! in progress) is refused, naming them. `--check` fetches the pinned commit into a scratch
//! repository only, so CI can tell that the pin exists. `goldens-pin` and `goldens-bless-pr`:
//! [`pin`], [`bless_pr`].

use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The data repository; `$SOLIDCORE_GOLDENS_URL` overrides it (a mirror, a local clone).
const URL: &str = "https://github.com/tangency-labs/solidcore-goldens.git";

fn url() -> String {
    env::var("SOLIDCORE_GOLDENS_URL").unwrap_or_else(|_| URL.to_string())
}

fn git(dir: &Path, args: &[&str]) -> io::Result<String> {
    let out = Command::new("git").arg("-C").arg(dir).args(args).output()?;
    if out.status.success() {
        Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
    } else {
        Err(io::Error::other(format!(
            "git {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr).trim()
        )))
    }
}

/// Fetches `sha` alone (depth 1) into the repository at `dir`.
fn fetch_sha(dir: &Path, sha: &str) -> io::Result<()> {
    git(dir, &["fetch", "--quiet", "--depth", "1", &url(), sha]).map(|_| ())
}

fn data_dir() -> PathBuf {
    match env::var_os(goldens_data::ENV) {
        Some(d) if !d.is_empty() => PathBuf::from(d),
        _ => goldens_data::repo_root().join(goldens_data::DEFAULT_DIR),
    }
}

/// `goldens-fetch` (`check`: `--check`).
pub fn fetch(check: bool) -> io::Result<()> {
    let sha = goldens_data::pin().map_err(io::Error::other)?;
    if check {
        let scratch = env::temp_dir().join(format!("goldens-check-{}", std::process::id()));
        fs::create_dir_all(&scratch)?;
        let result = git(&scratch, &["init", "--quiet", "--bare"])
            .and_then(|_| fetch_sha(&scratch, &sha))
            .and_then(|()| git(&scratch, &["rev-parse", "FETCH_HEAD"]));
        fs::remove_dir_all(&scratch)?;
        let got = result.map_err(|e| {
            io::Error::other(format!("the pinned sha {sha} is not on {}: {e}", url()))
        })?;
        if got != sha {
            return Err(io::Error::other(format!("fetched {got}, pinned {sha}")));
        }
        println!("goldens-fetch --check: {sha} is on {}", url());
        return Ok(());
    }
    let dir = data_dir();
    if dir.join(".git").exists() {
        let changed = git(&dir, &["status", "--porcelain"])?;
        if !changed.is_empty() {
            return Err(io::Error::other(format!(
                "{} has changes (a bless in progress?); commit them to solidcore-goldens or \
                 discard them first:\n{changed}",
                dir.display()
            )));
        }
        if git(&dir, &["rev-parse", "HEAD"]).ok().as_deref() == Some(sha.as_str()) {
            println!("goldens-fetch: {} is at {sha}", dir.display());
            return Ok(());
        }
    } else if dir.exists() && fs::read_dir(&dir)?.next().is_some() {
        return Err(io::Error::other(format!(
            "{} exists and is not a git checkout; move it aside",
            dir.display()
        )));
    } else {
        fs::create_dir_all(&dir)?;
        git(&dir, &["init", "--quiet"])?;
    }
    fetch_sha(&dir, &sha)?;
    git(
        &dir,
        &[
            "-c",
            "advice.detachedHead=false",
            "checkout",
            "--quiet",
            "--detach",
            "FETCH_HEAD",
        ],
    )?;
    let at = git(&dir, &["rev-parse", "HEAD"])?;
    if at != sha {
        return Err(io::Error::other(format!(
            "{} is at {at}, pinned {sha}",
            dir.display()
        )));
    }
    goldens_data::root().map_err(io::Error::other)?;
    println!("goldens-fetch: {} is at {sha}", dir.display());
    Ok(())
}

/// `goldens-pin <sha>`: writes `sha` (a full sha on solidcore-goldens) to `goldens.pin`,
/// after checking it is there, and fetches it.
pub fn pin(sha: &str) -> io::Result<()> {
    let sha = sha.trim().to_ascii_lowercase();
    if sha.len() != 40 || !sha.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(io::Error::other(format!(
            "goldens-pin: not a full sha: `{sha}`"
        )));
    }
    let path = goldens_data::repo_root().join(goldens_data::PIN_FILE);
    let old = fs::read_to_string(&path)?;
    fs::write(&path, format!("{sha}\n"))?;
    if let Err(e) = fetch(true).and_then(|()| fetch(false)) {
        fs::write(&path, old)?;
        return Err(e);
    }
    println!("goldens-pin: {} pins {sha}", path.display());
    Ok(())
}

/// `goldens-bless-pr <branch> <title> <body-file>`: the blessed changes in the data checkout
/// as a pull request on solidcore-goldens: a branch `<branch>` from the pin, one commit of
/// every change (the title its subject, the body file its message), pushed to the data
/// repository, and the pull request opened with `gh` against its `main`, the body file its
/// description (the justification the bless flow asks for). It prints the commit's sha;
/// pin it with `goldens-pin` once the pull request merges (its merge sha).
pub fn bless_pr(branch: &str, title: &str, body_file: &Path) -> io::Result<()> {
    let dir = data_dir();
    if !dir.join(".git").exists() {
        return Err(io::Error::other(format!(
            "{} is not a checkout; run `{}` and bless first",
            dir.display(),
            goldens_data::FETCH
        )));
    }
    let body = fs::read_to_string(body_file)?;
    if body.trim().is_empty() {
        return Err(io::Error::other(
            "the body file is empty: it carries the justification (cases, old and new values, why)",
        ));
    }
    let sha = goldens_data::pin().map_err(io::Error::other)?;
    if git(&dir, &["rev-parse", "HEAD"])? != sha {
        return Err(io::Error::other(format!(
            "{} is not at the pin {sha}; bless on the pinned data",
            dir.display()
        )));
    }
    let changed = git(&dir, &["status", "--porcelain"])?;
    if changed.is_empty() {
        return Err(io::Error::other(
            "nothing to propose: no changes in the data checkout",
        ));
    }
    println!("goldens-bless-pr: proposing\n{changed}");
    git(&dir, &["switch", "--quiet", "--create", branch])?;
    git(&dir, &["add", "--all", "--", "expected", "fixtures"])?;
    let message = format!("{title}\n\n{body}");
    git(&dir, &["commit", "--quiet", "--message", &message])?;
    let commit = git(&dir, &["rev-parse", "HEAD"])?;
    git(
        &dir,
        &[
            "push",
            "--quiet",
            &url(),
            &format!("HEAD:refs/heads/{branch}"),
        ],
    )?;
    let out = Command::new("gh")
        .args([
            "pr",
            "create",
            "--repo",
            "tangency-labs/solidcore-goldens",
            "--base",
            "main",
        ])
        .args(["--head", branch, "--title", title, "--body-file"])
        .arg(body_file)
        .output()?;
    if !out.status.success() {
        return Err(io::Error::other(format!(
            "gh pr create: {} (the branch {branch} is pushed at {commit})",
            String::from_utf8_lossy(&out.stderr).trim()
        )));
    }
    println!(
        "goldens-bless-pr: {}commit {commit} on {branch}; once it merges, \
         `cargo xtask goldens-pin <merge sha>` in your solidcore pull request",
        String::from_utf8_lossy(&out.stdout)
    );
    Ok(())
}
