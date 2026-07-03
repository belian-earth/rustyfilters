use extendr_api::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

static NUM_THREADS: AtomicUsize = AtomicUsize::new(1);
static POOL: Mutex<Option<rayon::ThreadPool>> = Mutex::new(None);

pub(crate) fn get_num_threads() -> usize {
    NUM_THREADS.load(Ordering::Relaxed)
}

/// `n == 0` means auto-detect: use all available cores.
fn set_num_threads(n: usize) {
    let n = if n == 0 {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(1)
    } else {
        n
    };
    NUM_THREADS.store(n, Ordering::Relaxed);
    if n > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .expect("failed to build thread pool");
        *POOL.lock().unwrap() = Some(pool);
    } else {
        *POOL.lock().unwrap() = None;
    }
}

/// Run closure on the rayon pool if threads > 1, otherwise run directly.
pub(crate) fn maybe_par<F, R>(f: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    let guard = POOL.lock().unwrap();
    match guard.as_ref() {
        Some(pool) => pool.install(f),
        None => f(),
    }
}

/// @noRd
/// @keywords internal
#[extendr]
fn rf_set_threads_rs(n: i32) {
    set_num_threads(n as usize);
}

/// @noRd
/// @keywords internal
#[extendr]
fn rf_get_threads_rs() -> i32 {
    get_num_threads() as i32
}

extendr_module! {
    mod threading;
    fn rf_set_threads_rs;
    fn rf_get_threads_rs;
}
