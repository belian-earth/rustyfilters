use extendr_api::prelude::*;

mod engine;
mod focal;
mod gaussian;
mod threading;

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyfilters;
    use focal;
    use gaussian;
    use threading;
}
