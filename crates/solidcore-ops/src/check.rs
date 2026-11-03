//! The kernel's validity check: `SmBrep::ValidatePointers`, then `AssertValid` over the whole
//! topology graph, reported a line per failed test.

use std::fmt::Write as _;

use solidcore::SolidError;
use solidcore::geom::assert_array::{AssertArray, AssertReport, obj_type_string};
use solidcore::topology::NamingScheme;
use solidcore::topology::brep::validation::{CheckLevel as KernelCheckLevel, check_brep};

use crate::body::Body;
use crate::session::with_body;
use crate::{Error, Result};


/// How many of solidcore's tests to run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckLevel {
    /// The fewest: topology and pointers, no geometry.
    Basic = 0,
    /// Adds the geometry checks (gaps between vertices, edges and faces against tolerances).
    Geometry = 1,
    /// Every test, the slow ones too.
    All = 2,
}

/// What solidcore's check found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CheckReport {
    /// Failed tests; 0 means solidcore finds the body valid.
    pub failures: usize,
    /// One line per failed test, such as
    /// `SmBrep::Test #82 Tolerance: Vertex tolerance relative to brep tolerance (offset=2, tolerance=0.01) [owner SmVertex]`.
    pub text: String,
}

impl CheckReport {
    /// Whether no test failed.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.failures == 0
    }

    /// The number of each failed test, line by line (`SmBrep::Test #82 …` is 82); a line without
    /// one is left out.
    #[must_use]
    pub fn failed_tests(&self) -> Vec<u32> {
        self.text
            .lines()
            .filter_map(|l| {
                let rest = l.split_once("Test #")?.1;
                let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
                digits.parse().ok()
            })
            .collect()
    }
}

impl std::fmt::Display for CheckReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_valid() { f.write_str("valid") } else { write!(f, "{} failed tests:\n{}", self.failures, self.text) }
    }
}

/// SMLib's `SM_SUCCESS`, the status `SmBrep::ValidatePointers` returns when it finds nothing.
const SM_SUCCESS: i32 = 1000;

/// Why [`Body::check`] could not run.
enum Refusal {
    /// A test the kernel has not ported (its [`SolidError::Unported`]).
    Unported(&'static str, &'static str),
    /// A body whose geometry owners' identities the check cannot know.
    Naming(NamingScheme),
}

impl crate::session::Answer for std::result::Result<solidcore::topology::brep::validation::BrepCheck, Refusal> {}

impl Body {
    /// Runs the kernel's validity check at `level`. A body that fails it is still `Ok`, with the
    /// failures in the report; `Err` means the check could not run.
    ///
    /// # Errors
    /// "not supported yet" while a test it reaches is unported, or for a body named by its
    /// maker's keys; those of [`with_body`].
    pub fn check(&self, level: CheckLevel) -> Result<CheckReport> {
        let lv = match level {
            CheckLevel::Basic => KernelCheckLevel::Level0,
            CheckLevel::Geometry => KernelCheckLevel::Level1,
            CheckLevel::All => KernelCheckLevel::Level2,
        };
        let c = match with_body("check", self, |k, b| {
            let naming = b.naming();
            if naming == NamingScheme::Unknown {
                return Ok(Err(Refusal::Naming(naming)));
            }
            match check_brep(k, b, lv) {
                Err(SolidError::Unported { upstream, unit }) => Ok(Err(Refusal::Unported(upstream, unit))),
                r => r.map(Ok),
            }
        })? {
            Ok(c) => c,
            Err(Refusal::Unported(upstream, unit)) => {
                return Err(Error::Kernel(format!("check: {upstream}: not supported yet: solidcore has not ported it (unit {unit})")));
            }
            Err(Refusal::Naming(n)) => {
                return Err(Error::Kernel(format!(
                    "check: not supported yet: the body names its geometry owners by its maker's keys (NamingScheme::{n:?}), which the check cannot know"
                )));
            }
        };
        let mut text = String::new();
        let mut failures = 0usize;
        collect(&c.pointer_reports, &mut text, &mut failures);
        if c.pointers_status != SM_SUCCESS && failures == 0 {
            failures += 1;
            let _ = writeln!(text, "SmBrep::ValidatePointers failed with status {}", c.pointers_status);
        }
        if failures > 0 {
            // The capi does not run AssertValid after a failed pointer check. check_brep runs it when
            // the status failed without a report; that result is dropped here, as the capi never
            // computes it.
            text.push_str("AssertValid skipped: ValidatePointers failed\n");
        } else {
            let Some(ok) = c.assert_valid else {
                return Err(Error::Kernel("check: AssertValid did not run after a clean pointer check".into()));
            };
            collect(&c.reports, &mut text, &mut failures);
            if !ok && failures == 0 {
                failures += 1;
                text.push_str("SmBrep::AssertValid failed without a report\n");
            }
        }
        Ok(CheckReport { failures, text })
    }
}

/// The capi's `text()`: an empty string reads as `?`.
fn text(s: &str) -> &str {
    if s.is_empty() { "?" } else { s }
}

/// The capi's `line()`: `<reporting type>::Test #<n> <name>: <message> [owner <type>]`, the owner
/// tag only when the report names the owner's type.
fn report_line(r: &AssertReport) -> String {
    let msg = r.format_log_message();
    let mut s = format!("{}::Test #{} {}: {}", text(obj_type_string(r.reporting_type)), r.test_index, text(r.name), text(&msg));
    if let Some(owner) = r.owner_type_string.filter(|o| !o.is_empty()) {
        s.push_str(" [owner ");
        s.push_str(owner);
        s.push(']');
    }
    s
}

/// The capi's `collect()`: appends the failed reports of `list`, one line each, and counts them;
/// SMLib's warnings ([`AssertReport::is_warning`]: `SmBrep` test 21, "Suspicious Shell Count", a
/// body of several faces with one shell each, is one) are neither.
fn collect(list: &AssertArray, out: &mut String, failures: &mut usize) {
    for r in list.reports().iter().filter(|r| !r.ok && !r.is_warning()) {
        *failures += 1;
        out.push_str(&report_line(r));
        out.push('\n');
    }
}
