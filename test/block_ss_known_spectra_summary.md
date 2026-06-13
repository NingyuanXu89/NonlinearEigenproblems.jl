# Block SS Known-Spectrum Benchmark Audit

## Purpose

This summary records the constructor conventions, contour solver interfaces,
existing coverage overlap, planned benchmark cases, stress-case expectations,
and narrow validation commands to use before adding reusable known-spectrum
tests for `contour_block_SS`.

## Constructor And Helper Conventions

- Use `PEP` for monomial polynomial NEPs.
- Use `DEP` for delay NEPs of the form `M(lambda) = -lambda*I + sum(A_i * exp(-tau_i*lambda))`.
- Use `SPMF_NEP` for sums of matrix coefficients and scalar/matrix functions.
- Use `Mder_NEP` and `Mder_Mlincomb_NEP` only for custom function-handle cases that do not fit `PEP`, `DEP`, or `SPMF_NEP`.
- Follow existing reusable example style: `nep_gallery` names and local helper functions are lowercase, such as `dep0`, `pep0`, `qep_fixed_eig`, `sine_nep`, and `diagonal_linear_pep`.

## Contour Solver Interfaces

The actual `contour_block_SS` call pattern uses the Unicode keyword `σ`, as
shown in existing tests and docs:

```julia
lambda, V = contour_block_SS(nep; σ=..., radius=..., N=..., k=..., K=...)
info = contour_block_SS_info(nep; σ=..., radius=..., N=..., k=..., K=...)
```

Typed and integrator forms are also available:

```julia
contour_block_SS(ComplexF64, nep; ...)
contour_block_SS(nep, MatrixTrapezoidal; ...)
contour_block_SS_info(ComplexF64, nep; ...)
contour_block_SS_info(nep, MatrixTrapezoidal; ...)
```

`contour_block_SS_info` returns diagnostics with these fields:
`lambda`, `V`, `singular_values`, `estimated_rank`, `rank_drop_tol`,
`capacity`, `number_returned`, `residuals`, and `inside_contour`.

## Existing Coverage

- `test/contour_block_SS.jl` covers circle and ellipse contours, `Shat_mode=:JSIAM`, residual smoke checks, and consistency between `contour_block_SS` and `contour_block_SS_info`.
- `test/beyn.jl` provides analogous Beyn contour and diagnostics coverage.
- `test/sim.jl` defines `diagonal_linear_pep(lambda_values)` and covers SIM region helpers, screening, probe validation, and contour decision policy.
- No existing test file provides a reusable known-spectrum benchmark suite for `contour_block_SS`.

## Benchmark Case Disposition

Implemented benchmark families:

- diagonal polynomial NEP;
- diagonal trigonometric-polynomial NEP;
- shifted diagonal exponential NEP;
- diagonal delay-type NEP with known Lambert-W roots;
- similarity-transformed diagonal NEP;
- ill-scaled diagonal NEP;
- defective or Jordan-like NEP;
- clustered eigenvalues near a contour boundary.

Planned benchmark families: none in the current known-spectrum split.

The delay benchmark uses `DEP` with `A = Diagonal(roots .* exp.(τ .* roots))`,
so each supplied root `r` satisfies
`-r + (r * exp(τ*r)) * exp(-τ*r) = 0`. Equivalently, for
`a_i = r_i * exp(τ*r_i)`, `r_i = W_0(τ*a_i) / τ` when it lies on the
principal Lambert-W branch. These roots are precomputed/chosen verified roots
for the constructed delay equation; the tests do not call
`SpecialFunctions.lambertw`.

Skipped cases: none.

The ill-scaled diagonal benchmark is diagnostics-oriented: it verifies block SS
diagnostic consistency, finite residuals, and known inside roots with looser
tolerances instead of requiring stronger behavior than the current solver
reports for deliberately varied component scales.

The Jordan-like and clustered-boundary benchmarks are stress cases: they verify
diagnostic visibility and documented sensitivity rather than requiring perfect
known-spectrum recovery.

Known overlap:

- The diagonal polynomial case partially overlaps with SIM's
  `diagonal_linear_pep(lambda_values)`, but that helper is scoped to SIM
  screening. The dedicated `diagonal_polynomial_pep` benchmark now covers
  block SS known-spectrum recovery directly.
- The remaining planned cases have no sufficient dedicated block SS
  known-spectrum coverage from the audit.

## Stress-Case Sensitivity

Clustered-boundary, defective or Jordan-like, and ill-scaled cases should not
force perfect recovery. These tests should verify that available diagnostics
expose sensitivity where possible, for example through residuals, rank-gap
behavior, capacity pressure, or `inside_contour` flags.

## Per-Problem `contour_block_SS` Results

Validated on 2026-06-13 by running `contour_block_SS` directly on each
known-spectrum problem with its configured contour and solver keywords.

| Problem | Constructor style | Result |
| --- | --- | --- |
| `diagonal_polynomial_pep` | `PEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `2.2e-16`. |
| `diagonal_trig_polynomial_nep` | `SPMF_NEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `2.6e-16`. |
| `shifted_diagonal_exponential_nep` | `SPMF_NEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `2.6e-16`. |
| `diagonal_delay_lambertw_nep` | `DEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `2.8e-16`. |
| `similarity_diagonal_nep` | `PEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `2.6e-16`. |
| `ill_scaled_diagonal_nep` | `PEP` | Returned 2 eigenvalues; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `1.6e-14`. |
| `jordan_like_nep` | `PEP` | Returned 3 eigenvalues; recovered all 3 desired inside roots, including the repeated defective target; accepted 0 listed outside roots; max residual about `8.1e-16`. |
| `clustered_boundary_nep` | `PEP` | Returned 4 eigenvalues, with 2 inside the contour; recovered all 2 desired inside roots; accepted 0 listed outside roots; max residual about `1.8e-9`. |

## Narrow Test Commands

Run only new or narrow relevant tests, not the full package test suite.

Preferred command for the later executable benchmark tests:

```bash
julia --project=. test/block_ss_known_spectra.jl
```

This repository also has `test/Project.toml`, so use the test project only if
the root-project command cannot run the narrow test structurally:

```bash
julia --project=test test/block_ss_known_spectra.jl
julia --project=test test/runtests.jl block_ss_known_spectra
```

Only run `Pkg.test()` if narrow commands cannot run structurally, and explain
why that broader command was necessary.
