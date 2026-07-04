use extendr_api::prelude::*;

mod bilateral;
mod convolve;
mod engine;
mod focal;
mod gaussian;
mod guided;
mod speckle;
mod threading;

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyfilters;
    use bilateral;
    use convolve;
    use focal;
    use gaussian;
    use guided;
    use speckle;
    use threading;
}
