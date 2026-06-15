# Block SS known-spectrum benchmarks

The test suite includes a known-spectrum benchmark set for
[`contour_block_SS`](@ref) and [`contour_block_SS_info`](@ref). These tests are
intended to exercise block SS recovery and diagnostics on problems where the
inside and outside eigenvalues are known in advance.

The executable tests live in `test/block_ss_known_spectra.jl`, with reusable
problem constructors in `test/nep_benchmarks/block_ss_known_spectra_problems.jl`.
The longer maintainer audit is kept in
`test/block_ss_known_spectra_summary.md`.

## Benchmark families

The current suite covers:

- diagonal polynomial NEPs;
- diagonal trigonometric-polynomial NEPs;
- shifted diagonal exponential NEPs;
- delay/Lambert-W style NEPs;
- similarity-transformed diagonal NEPs;
- ill-scaled diagonal NEPs;
- defective or Jordan-like NEPs;
- clustered eigenvalues near a contour boundary.

The clean cases validate recovery of the known inside eigenvalues and rejection
of listed outside values. The stress cases are deliberately more sensitive:
they check that block SS diagnostics such as residuals, estimated rank,
capacity, and inside-contour flags remain informative rather than requiring
perfect recovery under ill scaling, defective structure, or boundary clustering.

## Narrow validation command

Run this suite directly when changing block SS behavior or diagnostics:

```bash
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' block_ss_known_spectra.jl
```

This is a targeted benchmark-style validation, not a replacement for the
smaller smoke tests in `test/contour_block_SS.jl`.
