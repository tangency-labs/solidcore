//! The session: one kernel for the process (SMLib's context is process-global too, and some of
//! its state, such as the fillet solver's latched step, lives for the whole process), every call
//! under its lock and a panic guard.
//!
//! - **Fixes:** the kernel runs with the opt-in fixes the application names ([`set_fixes`], by
//!   `UpstreamFixes` flag name); which ones is the application's choice, never ops'.
//! - **Panics:** a panic inside a call is caught; the kernel, which it may have left half-updated,
//!   is replaced by a fresh one with the same fixes and the call fails with [`Error::Internal`].
//! - **Failures:** an application may install a callback ([`set_failure_hook`]) that is handed a
//!   failed call's name, why it failed and its input bodies as binary `.smb` bytes. Ops itself
//!   writes nothing and logs nothing.
//! - **The call running now:** [`current_call`] reads it without a lock, for crash reports.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};

use solidcore::kernel::{Kernel, KernelOptions, UpstreamFixes};
use solidcore::{Brep, SolidError};

use crate::body::Body;
use crate::{Error, Result};

static KERNEL: Mutex<Option<Kernel>> = Mutex::new(None);
static FIXES: Mutex<Option<UpstreamFixes>> = Mutex::new(None);

/// A failed call, as the failure callback sees it.
#[derive(Debug)]
pub struct Failure<'a> {
    /// The call's name.
    pub call: &'a str,
    /// Why it failed: its error, or the refusal its answer carries.
    pub why: &'a str,
    /// Its input bodies by the names the call gives them, each as a binary `.smb` (or why it
    /// could not be written), as they are now: a call that edits a body in place leaves it edited.
    pub inputs: Vec<(&'static str, std::result::Result<Vec<u8>, String>)>,
}

type Hook = Box<dyn Fn(&Failure<'_>) + Send + Sync>;
static HOOK: Mutex<Option<Hook>> = Mutex::new(None);

/// A call's input bodies, by the names its failure report gives them, in the order it takes them.
pub type Inputs<'a> = &'a [(&'static str, &'a Body)];

/// A flag name `UpstreamFixes` does not have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnknownFix(pub String);

impl std::fmt::Display for UnknownFix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "solidcore has no upstream fix named {}", self.0)
    }
}

impl std::error::Error for UnknownFix {}

/// The kernel runs with exactly the fixes named in `on` (each an `UpstreamFixes` field name;
/// every other flag off) from the next call on. Nothing changes if a name is unknown.
///
/// # Errors
/// [`UnknownFix`] naming the first name `UpstreamFixes` does not have.
pub fn set_fixes<'a>(on: impl IntoIterator<Item = &'a str>) -> std::result::Result<(), UnknownFix> {
    let mut fixes = UpstreamFixes::default();
    for name in on {
        *fixes.flag_mut(name).ok_or_else(|| UnknownFix(name.to_owned()))? = true;
    }
    *FIXES.lock().unwrap_or_else(PoisonError::into_inner) = Some(fixes);
    let mut kernel = KERNEL.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some(k) = kernel.as_mut() {
        k.fixes = fixes;
    }
    Ok(())
}

/// The names of the fixes the kernel runs with now.
#[must_use]
pub fn fixes() -> Vec<&'static str> {
    current_fixes().flags().filter(|(_, on)| *on).map(|(name, _)| name).collect()
}

fn current_fixes() -> UpstreamFixes {
    FIXES.lock().unwrap_or_else(PoisonError::into_inner).unwrap_or_default()
}

/// Installs the failure callback (`None` removes it). It runs under the kernel's lock and
/// must not call ops.
pub fn set_failure_hook(hook: Option<Hook>) {
    *HOOK.lock().unwrap_or_else(PoisonError::into_inner) = hook;
}

fn fresh() -> Kernel {
    let mut k = Kernel::new(&KernelOptions::default());
    k.fixes = current_fixes();
    k
}

thread_local! {
    /// How many kernel calls this thread is inside: a panic now is caught ([`run_on`]).
    static IN_CALL: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

/// Is this thread inside a kernel call, whose panic the session catches? (Crash reports mark
/// such a panic recovered.)
#[must_use]
pub fn in_call() -> bool {
    IN_CALL.with(std::cell::Cell::get) > 0
}

/// The longest call name [`current_call`] gives.
const NAME_MAX: usize = 96;
static NAME: [AtomicU8; NAME_MAX] = [const { AtomicU8::new(0) }; NAME_MAX];
static NAME_LEN: AtomicUsize = AtomicUsize::new(0);

/// `call` is running now (`None`: none is). One at a time, under the kernel's lock.
fn set_current(call: Option<&str>) {
    NAME_LEN.store(0, Ordering::SeqCst);
    let Some(call) = call else { return };
    let bytes = call.as_bytes();
    let n = bytes.len().min(NAME_MAX);
    for (slot, b) in NAME.iter().zip(bytes) {
        slot.store(*b, Ordering::Relaxed);
    }
    NAME_LEN.store(n, Ordering::SeqCst);
}

/// The kernel call running now, if any. Takes no lock, for crash reports, so it may be cut short
/// or mixed with the next call's name if read just as one starts.
#[must_use]
pub fn current_call() -> Option<String> {
    let n = NAME_LEN.load(Ordering::SeqCst);
    if n == 0 {
        return None;
    }
    let bytes: Vec<u8> = NAME.iter().take(n).map(|b| b.load(Ordering::Relaxed)).collect();
    Some(String::from_utf8_lossy(&bytes).into_owned())
}

/// How a call's value reads: whether it is a refusal, for the calls whose answer can carry one.
pub trait Answer {
    /// The refusal this value carries, if it is one.
    fn refusal(&self) -> Option<String> {
        None
    }
}

impl<T> Answer for Result<T> {
    fn refusal(&self) -> Option<String> {
        self.as_ref().err().map(Error::to_string)
    }
}

/// The values that are no refusal.
macro_rules! plain {
    ($($t:ty),* $(,)?) => { $(impl Answer for $t {})* };
}
plain!((), bool, f64, usize, [usize; 3], Vec<u8>, Brep, crate::check::CheckReport);

impl<T> Answer for Option<T> {}

/// Runs `f` on the process's kernel ([`run_on`] with no inputs).
///
/// # Errors
/// As [`run_on`].
pub fn run<T: Answer>(call: &'static str, f: impl FnOnce(&mut Kernel) -> solidcore::Result<T>) -> Result<T> {
    run_on(call, &[], f)
}

/// Runs `f` on the process's kernel, for a call that reads the bodies `inputs`. Should it fail
/// while a failure callback is installed, the callback gets them.
///
/// # Errors
/// [`Error::Kernel`] for the kernel's error, [`Error::Internal`] for a panic.
pub fn run_on<T: Answer>(call: &'static str, inputs: Inputs<'_>, f: impl FnOnce(&mut Kernel) -> solidcore::Result<T>) -> Result<T> {
    struct Depth;
    impl Drop for Depth {
        fn drop(&mut self) {
            IN_CALL.with(|d| d.set(d.get().saturating_sub(1)));
        }
    }
    let mut guard = KERNEL.lock().unwrap_or_else(PoisonError::into_inner);
    let k = guard.get_or_insert_with(fresh);
    set_current(Some(call));
    IN_CALL.with(|d| d.set(d.get() + 1));
    let depth = Depth;
    let caught = catch_unwind(AssertUnwindSafe(|| f(k)));
    drop(depth);
    let r = match caught {
        Ok(Ok(v)) => Ok(v),
        Ok(Err(e)) => Err(Error::Kernel(format!("{call}: {e}"))),
        Err(p) => {
            *guard = Some(fresh());
            let why = p.downcast_ref::<&str>().map(|s| (*s).to_owned()).or_else(|| p.downcast_ref::<String>().cloned()).unwrap_or_default();
            Err(Error::Internal(format!("{call}: {why}")))
        }
    };
    set_current(None);
    let why = match &r {
        Ok(v) => v.refusal(),
        Err(e) => Some(e.to_string()),
    };
    if let Some(why) = why {
        // The kernel the call ran on, or the fresh one that replaced it after a panic.
        failed(guard.get_or_insert_with(fresh), call, &why, inputs);
    }
    r
}

/// Hands a failed call to the failure callback, if one is installed; its inputs are written only
/// then.
fn failed(k: &Kernel, call: &str, why: &str, inputs: Inputs<'_>) {
    let hook = HOOK.lock().unwrap_or_else(PoisonError::into_inner);
    let Some(hook) = hook.as_ref() else { return };
    let inputs = inputs
        .iter()
        .map(|&(what, body)| (what, crate::file::bytes(k, &body.lock(), false).map_err(|e| e.to_string())))
        .collect();
    hook(&Failure { call, why, inputs });
}

/// A call the kernel answers with a status, any status but success standing for no answer: an
/// error that is exactly a kernel status ([`SolidError::sm_status`]) is `None`. Every other
/// error (unported code, invalid input to an adapter, a missing record...) fails the call
/// instead of hiding as a skipped entity.
///
/// # Errors
/// The kernel's error when it stands for no status.
pub fn answer<T>(r: solidcore::Result<T>) -> solidcore::Result<Option<T>> {
    match r {
        Ok(v) => Ok(Some(v)),
        Err(e) if e.sm_status().is_some() => Ok(None),
        Err(e) => Err(e),
    }
}

/// Runs `f` on the kernel and `body`, locking the kernel first.
///
/// # Errors
/// [`Error::Internal`] for a body an earlier panic left unfinished; those of [`run_on`].
pub fn with_body<T: Answer>(call: &'static str, body: &Body, f: impl FnOnce(&mut Kernel, &mut Brep) -> solidcore::Result<T>) -> Result<T> {
    // A body a panic interrupted may be half-edited: it is never used again.
    if body.is_poisoned() {
        return Err(Error::Internal(format!("{call}: the body was left unfinished by an earlier panic")));
    }
    run_on(call, &[("brep", body)], |k| f(k, &mut body.lock()))
}

/// A deliberate panic inside a kernel call, which the session catches as it catches any
/// ([`Error::Internal`], the kernel replaced by a fresh one): for an application's crash-report
/// check.
///
/// # Errors
/// Always [`Error::Internal`].
#[allow(clippy::panic)]
pub fn crash_test() -> Result<()> {
    run("crash test", |_| -> std::result::Result<(), SolidError> { std::panic::panic_any(String::from("crash test: an internal error in the kernel")) })
}
