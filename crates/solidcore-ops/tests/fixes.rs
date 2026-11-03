//! The fix set, taken by name: its own test binary, as the fix set is the process kernel's.

use solidcore_ops::session;

/// The fixes are the ones named, by `UpstreamFixes` field name; an unknown name changes nothing.
#[test]
fn fixes_are_taken_by_name() {
    let name = solidcore::kernel::UpstreamFixes::NAMES[0];
    let before = session::fixes();
    assert!(session::set_fixes(["no_such_fix"]).is_err());
    assert_eq!(session::fixes(), before);
    session::set_fixes([name]).unwrap();
    assert_eq!(session::fixes(), vec![name]);
    assert!(session::run("fixes", |k| Ok(k.fixes.flags().any(|(n, on)| n == name && on))).unwrap());
    session::set_fixes([]).unwrap();
    assert!(session::fixes().is_empty());
}
