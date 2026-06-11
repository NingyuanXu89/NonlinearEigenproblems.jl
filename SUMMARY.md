# SIM Workflow Change Summary

This branch adds a modular Spectral Indicator Method (SIM) workflow for
screening regions before applying existing contour eigensolvers. The original
repository had Beyn and block SS contour solvers, plus shared contour
quadrature helpers, but it did not have region abstractions, SIM screening,
enclosing-contour conversion, contour diagnostic return types, or a policy
layer connecting screening with solver diagnostics.

## Architecture

The implementation keeps the new workflow layered:

1. **Region and enclosing-contour helpers**
   - Added in `src/method_sim.jl`.
   - New public types:
     - `SquareRegion`
     - `RectangularRegion`
     - `EllipseContour`
   - New public helpers:
     - `enclosing_contour`
     - `contour_parameters`
     - `inside_region`
     - `subdivide`
   - These helpers are independent of the contour solvers and do not modify
     `src/method_contour_common.jl`.

2. **SIM screening**
   - Added in `src/method_sim.jl`.
   - New public functions:
     - `sim_indicator`
     - `sim_screen`
     - `sim_screen_regions`
   - SIM computes the Xi-Sun projection-ratio indicator for active/inactive
     region screening.
   - SIM does not solve eigenvalue problems and does not count eigenvalues.
   - Thresholding is applied only by `sim_screen`; `sim_indicator` returns the
     raw screening scalar.

3. **Contour diagnostics**
   - Added non-breaking info variants:
     - `contour_beyn_info`
     - `contour_block_SS_info`
   - Existing public solver calls still return `(lambda, V)`:
     - `contour_beyn`
     - `contour_block_SS`
   - Beyn diagnostics report singular values of the `A0` moment matrix.
   - Block SS diagnostics report singular values of the block Hankel moment
     matrix.
   - Rank/capacity diagnostics are separate from SIM screening.

4. **Policy helper**
   - Added `SIMContourDecision` and `sim_contour_decision`.
   - The decision helper combines a `SIMResult` with optional contour
     diagnostics.
   - It is only a policy layer; it does not run SIM or contour solvers.

5. **Demo workflow**
   - Added `docs/src/sim_workflow_demo.jl`.
   - The demo uses a deterministic diagonal toy `PEP` with known eigenvalues.
   - It shows:
     - region definition,
     - subdivision,
     - SIM screening,
     - enclosing-contour conversion,
     - `contour_beyn_info`,
     - `inside_region` filtering,
     - policy decisions,
     - optional Newton refinement.
   - Added `test/sim_workflow_demo.jl` so the demo remains executable.

## Public API Additions

The following names are exported through `NEPSolver` and re-exported by the
top-level package:

- `SquareRegion`
- `RectangularRegion`
- `EllipseContour`
- `SIMResult`
- `SIMContourDecision`
- `enclosing_contour`
- `contour_parameters`
- `inside_region`
- `subdivide`
- `sim_indicator`
- `sim_screen`
- `sim_screen_regions`
- `sim_contour_decision`
- `contour_beyn_info`
- `contour_block_SS_info`
- `ContourBeynInfo`
- `ContourBlockSSInfo`

The existing `contour_beyn` and `contour_block_SS` APIs are preserved.

## Numerical Notes

- Rectangular and square regions are axis-aligned in the complex plane.
- `enclosing_contour(region; shape=:circle)` returns the circumscribed circle,
  with radius `sqrt(half_width^2 + half_height^2)`.
- `enclosing_contour(region; shape=:ellipse)` returns an axis-aligned ellipse
  with the region half-width and half-height as radii.
- SIM screening currently uses the circumscribed circle. For rectangular
  regions, this can intentionally produce false positives outside the
  rectangle; callers should filter final candidates with `inside_region`.
- The SIM indicator uses the Xi-Sun half-grid projection-ratio criterion:
  `norm(Pf_N ./ Pf_(N/2)) / sqrt(n)`.
- Division safeguards treat zero-over-zero projection components as zero and
  nonzero-over-zero components as `Inf`.

## Tests and Documentation

Added or extended tests:

- `test/sim.jl`
  - region construction and validation,
  - enclosing contours,
  - subdivision,
  - SIM active/inactive screening,
  - deterministic probe behavior,
  - decision-policy behavior.
- `test/beyn.jl`
  - `contour_beyn_info` compatibility and rank diagnostics.
- `test/contour_block_SS.jl`
  - `contour_block_SS_info` compatibility and rank diagnostics.
- `test/sim_workflow_demo.jl`
  - executable coverage for the demo script.

Documentation updates:

- `docs/src/methods.md`
  - includes the new public API in the contour-method docs list.
- `docs/src/tutorial_contour.md`
  - adds a short SIM screening workflow section.
- `docs/src/sim_workflow_demo.jl`
  - provides a deterministic runnable example.

## Verification

Focused verification used during development:

```julia
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' 'sim|beyn|contour_block_SS'
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("docs/src/sim_workflow_demo.jl")'
```

`git diff --check` was also run to catch whitespace issues.

Known unrelated baseline failures mentioned before this work, including
blocknewton, cd player native, and Jacobi-Davidson/Effenberger, were not
addressed.
