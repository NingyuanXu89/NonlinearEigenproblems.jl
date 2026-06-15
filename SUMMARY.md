# SIM Workflow Change Summary

This branch adds a region-based Spectral Indicator Method (SIM) workflow for
cheaply screening contours before running contour eigensolvers. The workflow now
supports the legacy Xi-Sun screen, a low-cost moment-norm screen for subdivision,
method-aware diagnostics and decisions, local seeded contour probes, and an
auto-subdivision helper for preparing smaller block SS recovery regions.

## Architecture

1. **Region and plotting helpers**
   - Added `SquareRegion`, `RectangularRegion`, and `EllipseContour`.
   - Added `enclosing_contour`, `contour_parameters`, `inside_region`, and
     `subdivide`.
   - Added boundary collectors for plotting:
     `region_boundary`, `contour_boundary`, `collect_region_boundaries`, and
     `collect_contour_boundaries`.

2. **SIM screening**
   - Added `sim_indicator`, `sim_screen`, and `sim_screen_regions`.
   - `method=:xisun` preserves the legacy Xi-Sun / pmCIM projection-ratio
     screen.
   - `method=:moment_norm` computes low-order shifted contour moments with a
     fixed probe block and uses scaled moment norms as a screening heuristic.
   - `SIMResult` keeps the legacy fields and adds `diagnostics` for method,
     seed, moment strengths/norms/singular values, and contour metadata.

3. **Contour diagnostics and local seeds**
   - Added non-breaking info variants:
     `contour_beyn_info` and `contour_block_SS_info`.
   - Existing `contour_beyn` and `contour_block_SS` still return `(lambda, V)`.
   - `seed=10` is now exposed for Beyn and block SS random probe generation.
     Probes use local `MersenneTwister(seed)` objects and do not mutate Julia's
     global RNG state.

4. **Decision and subdivision policy**
   - Added `SIMContourDecision` and `sim_contour_decision`.
   - Inactive SIM decisions are treated as heuristic screens, not certificates
     of emptiness.
   - `method=:moment_norm` supports borderline handling and a
     `target_capacity` decision path for accepting or subdividing verified
     regions.
   - Added `sim_subdivide_active_regions`, which auto-subdivides an initial
     region with `method=:moment_norm`, verifies active regions with block SS,
     and returns accepted regions. With `return_trace=true`, it also returns
     plot-friendly path and decision metadata compatible with
     `collect_region_boundaries`.

5. **Moment-norm recovery workflow**
   - `docs/src/sim_workflow_demo.jl` now demonstrates the intended pipeline:
     auto-subdivide until each accepted region is roughly below `k*K`;
     enclose each accepted region by a circle or ellipse; run
     `contour_block_SS_info` with extra recovery moments
     `K' = K + additional_moments`; filter candidates with `inside_region`;
     refine candidates with Newton; de-duplicate with a user-chosen tolerance;
     and return solved eigenvalues/eigenvectors as `solved_λ` and `solved_V`.

## Public API Additions

The following names are exported through `NEPSolver` and re-exported by the
top-level package:

- `SquareRegion`
- `RectangularRegion`
- `EllipseContour`
- `SIMResult`
- `SIMContourDecision`
- `ContourBeynInfo`
- `ContourBlockSSInfo`
- `enclosing_contour`
- `contour_parameters`
- `inside_region`
- `subdivide`
- `region_boundary`
- `contour_boundary`
- `collect_region_boundaries`
- `collect_contour_boundaries`
- `sim_indicator`
- `sim_screen`
- `sim_screen_regions`
- `sim_contour_decision`
- `sim_subdivide_active_regions`
- `contour_beyn_info`
- `contour_block_SS_info`

The public return values of `contour_beyn` and `contour_block_SS` are preserved.

## Numerical Notes

- Rectangular and square regions are axis-aligned in the complex plane.
- `enclosing_contour(region; shape=:circle)` returns the circumscribed circle.
- `enclosing_contour(region; shape=:ellipse)` returns an axis-aligned ellipse
  with the region corners on the boundary.
- SIM screening uses enclosing contours and can produce false positives outside
  rectangular regions; final candidates should be filtered with `inside_region`.
- `method=:xisun` is a fast legacy projection-ratio screen based on the zeroth
  contour moment.
- `method=:moment_norm` is a low-order multi-moment screening heuristic. It does
  not recover eigenpairs, does not form block SS Hankel pencils, and does not
  prove that inactive regions are empty.
- Reliable eigenvalues should still be recovered with Beyn or block SS and
  validated with residuals, contour/subdivision changes, `N`, `k`, `K`, extra
  recovery moments, and optional Newton refinement.

## Tests and Documentation

Added or extended tests:

- `test/sim.jl`
  - region helpers and boundary collection,
  - legacy Xi-Sun screening,
  - moment-norm diagnostics and seed behavior,
  - method-aware decision policy,
  - auto-subdivision and trace metadata.
- `test/beyn.jl`
  - `contour_beyn_info` compatibility,
  - local seed reproducibility and global RNG isolation.
- `test/contour_block_SS.jl`
  - `contour_block_SS_info` compatibility,
  - local seed reproducibility and global RNG isolation.
- `test/sim_workflow_demo.jl`
  - executable coverage for the moment-norm workflow demo.

Documentation updates:

- `docs/src/methods.md`
  - includes the new public API in the contour-method docs list.
- `docs/src/tutorial_contour.md`
  - documents Xi-Sun and moment-norm SIM as heuristic screens, including
    auto-subdivision and trace metadata for plotting.
- `docs/src/block_ss_known_spectra.md`
  - summarizes the known-spectrum block SS benchmark families and narrow
    validation command.
- `docs/src/sim_workflow_demo.jl`
  - provides a deterministic runnable auto-subdivision and block SS recovery
    example.

## Verification

Focused verification used during development:

```julia
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' sim.jl
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' beyn.jl
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' contour_block_SS.jl
julia --project=test -e 'push!(LOAD_PATH, pwd()); include("test/runtests.jl")' sim_workflow_demo.jl
```

`git diff --check` was also run on edited files. Full `Pkg.test()` was not run
because the repository has unrelated long-running or known baseline failures.
