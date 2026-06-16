For each region in the BFS queue:

1. Run sim_screen(method=:moment_norm)

   This computes a moment-norm indicator for the enclosing circle of the region.
   The indicator is a heuristic measure of whether T(λ)^(-1) has eigenvalue
   singularities inside/near the contour.

2. Call sim_contour_decision(sim)

   A. inactive and far below threshold
      condition:
        indicator < threshold / borderline_factor

      meaning:
        moment strength is very small; the region is treated as empty/inactive.

      → action :skip
      → reason :moment_norm_inactive
      → status :skipped_inactive


   B. inactive but borderline
      condition:
        threshold / borderline_factor <= indicator <= threshold

      meaning:
        moment strength is weak but not negligible. This may be a false negative,
        for example from a small residue, cancellation, eigenvalue near the edge,
        or insufficient quadrature accuracy.

      → action :subdivide_or_rerun
      → reason :weak_moment_strength
      → status :subdivided_borderline


   C. active
      condition:
        indicator > threshold

      meaning:
        SIM detects nontrivial moment strength, so the region may contain one or
        more eigenvalues.

      → action :run_ss
      → reason :moment_detected


3. For active regions, run contour_block_SS_info

   This estimates eigenvalues/eigenvectors inside the enclosing contour and
   returns rank/capacity/residual diagnostics.

4. Call sim_contour_decision(sim, info; target_capacity)

   D. number_returned == 0
      meaning:
        SIM said the region was active, but block SS did not return eigenpairs.
        Possible causes include weak/ill-conditioned eigenvalues, bad rank
        estimation, poor quadrature, or an eigenvalue too close to the contour.

      → action :subdivide_or_rerun
      → reason :active_no_eigenpairs
      → status :subdivided_unreliable


   E. residuals too large, if residual_tol is used
      condition:
        any returned residual > residual_tol

      meaning:
        block SS returned candidate eigenpairs, but they do not satisfy
        T(λ)v ≈ 0 accurately enough. Possible causes include too small N,
        near-singular quadrature nodes, eigenvalues near the contour, bad
        conditioning, or too many modes inside one contour.

      → action :subdivide_or_rerun
      → reason :poor_residuals
      → status :subdivided_unreliable


   F. unclear rank gap
      condition:
        the singular-value drop used for rank estimation is not clean

      meaning:
        block SS cannot confidently determine the numerical rank / number of
        eigen-directions inside the contour. Possible causes include clustered
        eigenvalues, multiple eigenvalues, eigenvalues near the contour boundary,
        inappropriate rank_drop_tol, or insufficient quadrature accuracy.

      → action :subdivide_or_rerun
      → reason :unclear_rank_gap
      → status :subdivided_unreliable


   G. estimated_rank <= target_capacity
      meaning:
        the region appears to contain a manageable number of eigen-directions,
        and the rank estimate is reliable enough.

      → action :accept_region
      → reason :moment_rank_within_target
      → status :accepted


   H. estimated_rank > target_capacity
      meaning:
        the region appears to contain too many eigen-directions for the requested
        target capacity. Subdivision is used to split the spectrum into smaller
        pieces.

      → action :subdivide_or_rerun
      → reason :moment_rank_above_target
      → status :subdivided_over_target
