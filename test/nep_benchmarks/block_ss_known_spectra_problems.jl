# Helper definitions for known-spectrum block SS benchmark problems.
#
# This file is intended to be loaded with include(...). Later generator
# functions should return the named tuple created by
# block_ss_known_spectrum_problem. These helpers are test utilities only and do
# not change solver behavior.

using LinearAlgebra
using NonlinearEigenproblems

_as_named_tuple(nt::NamedTuple) = nt
_as_named_tuple(pairs) = (; pairs...)

_complex_vector(values) = collect(complex.(collect(values)))

"""
    block_ss_known_spectrum_problem(; name, nep, inside, outside, σ, radius,
                                     solver_kwargs, expected_behavior, notes)

Create the shared named-tuple representation used by known-spectrum block SS
benchmark generators.
"""
function block_ss_known_spectrum_problem(;
    name,
    nep,
    inside,
    outside=ComplexF64[],
    σ=0.0 + 0.0im,
    radius=1.0,
    solver_kwargs=(;),
    expected_behavior=:clean,
    notes=""
)
    return (;
        name=name,
        nep=nep,
        inside=_complex_vector(inside),
        outside=_complex_vector(outside),
        σ=σ,
        radius=radius,
        solver_kwargs=_as_named_tuple(solver_kwargs),
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

"""
    block_ss_solver_kwargs(problem; kwargs...)

Return contour solver keyword arguments for a known-spectrum problem. Explicit
keyword overrides take precedence over the problem's solver defaults.
"""
function block_ss_solver_kwargs(problem; kwargs...)
    return merge((; σ=problem.σ, radius=problem.radius),
                 problem.solver_kwargs,
                 (; kwargs...))
end

_radius_axes(radius::Number) = (radius, radius)
function _radius_axes(radius)
    length(radius) == 2 || error("radius must be a scalar or a length-two collection")
    return (radius[1], radius[2])
end

"""
    inside_ellipse(z; σ, radius)

Return true when `z` lies inside or on the ellipse used by the contour solvers.
"""
function inside_ellipse(z; σ, radius)
    radius1 = _radius_axes(radius)
    return (real(z - σ) / radius1[1])^2 + (imag(z - σ) / radius1[2])^2 <= 1
end

inside_problem_contour(problem, z) = inside_ellipse(z; σ=problem.σ, radius=problem.radius)

_match_tol(expected, atol, rtol) = atol + rtol * abs(expected)

"""
    match_eigenvalues(computed, expected; atol=1e-8, rtol=1e-8)

Greedily match expected eigenvalues to distinct computed eigenvalues.
"""
function match_eigenvalues(computed, expected; atol=1e-8, rtol=1e-8)
    computed_values = _complex_vector(computed)
    expected_values = _complex_vector(expected)

    distances = fill(Inf, length(expected_values))
    computed_indices = Vector{Union{Nothing,Int}}(nothing, length(expected_values))
    used = falses(length(computed_values))

    for (i, expected_value) in pairs(expected_values)
        best_index = nothing
        best_distance = Inf
        for (j, computed_value) in pairs(computed_values)
            used[j] && continue
            distance = abs(computed_value - expected_value)
            if distance < best_distance
                best_distance = distance
                best_index = j
            end
        end

        if best_index !== nothing
            distances[i] = best_distance
            computed_indices[i] = best_index
            used[best_index] = true
        end
    end

    matched = all(i -> distances[i] <= _match_tol(expected_values[i], atol, rtol),
                  eachindex(expected_values))
    return (; matched=matched, distances=distances, computed_indices=computed_indices)
end

function unmatched_expected(computed, expected; atol=1e-8, rtol=1e-8)
    expected_values = _complex_vector(expected)
    matches = match_eigenvalues(computed, expected_values; atol=atol, rtol=rtol)
    unmatched = Bool[
        !(matches.distances[i] <= _match_tol(expected_values[i], atol, rtol))
        for i in eachindex(expected_values)
    ]
    return expected_values[unmatched]
end

function accepted_outside_values(computed, outside; atol=1e-8, rtol=1e-8)
    outside_values = _complex_vector(outside)
    matches = match_eigenvalues(computed, outside_values; atol=atol, rtol=rtol)
    accepted = Bool[
        matches.distances[i] <= _match_tol(outside_values[i], atol, rtol)
        for i in eachindex(outside_values)
    ]
    return outside_values[accepted]
end

block_ss_eigenvalues(result::Tuple) = result[1]
function block_ss_eigenvalues(result)
    hasproperty(result, :lambda) || error("result does not contain eigenvalues")
    return getproperty(result, :lambda)
end

block_ss_vectors(result::Tuple) = result[2]
function block_ss_vectors(result)
    hasproperty(result, :V) || error("result does not contain eigenvectors")
    return getproperty(result, :V)
end

_optional_property(info, name) = hasproperty(info, name) ? getproperty(info, name) : nothing

function block_ss_diagnostics(info)
    return (;
        singular_values=_optional_property(info, :singular_values),
        estimated_rank=_optional_property(info, :estimated_rank),
        rank_drop_tol=_optional_property(info, :rank_drop_tol),
        capacity=_optional_property(info, :capacity),
        number_returned=_optional_property(info, :number_returned),
        residuals=_optional_property(info, :residuals),
        inside_contour=_optional_property(info, :inside_contour),
    )
end

const DIAGONAL_POLYNOMIAL_PEP_ROOTS =
    ComplexF64[-0.35 + 0.12im, 0.28 - 0.2im, 1.4 + 0.1im, -1.3 - 0.2im]

"""
    diagonal_polynomial_pep([roots]; kwargs...)

Return a diagonal linear `PEP` with known eigenvalues at `roots`.
"""
function diagonal_polynomial_pep(
    roots=DIAGONAL_POLYNOMIAL_PEP_ROOTS;
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:clean,
    notes="Diagonal linear polynomial PEP with two default roots inside the contour."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    nep = PEP([Matrix(Diagonal(-root_values)), Matrix{ComplexF64}(I, n, n)])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="diagonal_polynomial_pep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const DIAGONAL_TRIG_POLYNOMIAL_NEP_ROOTS =
    ComplexF64[-0.32 + 0.1im, 0.24 - 0.22im, 1.25 + 0.1im, -1.2 - 0.15im]

"""
    diagonal_trig_polynomial_nep([roots]; kwargs...)

Return a diagonal `SPMF_NEP` with entries
`m_j(λ) = (λ - roots[j]) + trig_scale * sin(λ - roots[j])`.
"""
function diagonal_trig_polynomial_nep(
    roots=DIAGONAL_TRIG_POLYNOMIAL_NEP_ROOTS;
    trig_scale=0.05,
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:clean,
    notes="Diagonal trigonometric-polynomial SPMF_NEP with roots away from the contour boundary."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    identity = Matrix{ComplexF64}(I, n, n)
    nep = SPMF_NEP(
        [
            identity,
            Matrix(Diagonal(-root_values)),
            Matrix(Diagonal(trig_scale .* cos.(root_values))),
            Matrix(Diagonal(-trig_scale .* sin.(root_values))),
        ],
        [S -> S, S -> one(S), S -> sin(S), S -> cos(S)]
    )

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="diagonal_trig_polynomial_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const SHIFTED_DIAGONAL_EXPONENTIAL_NEP_ROOTS =
    ComplexF64[-0.3 + 0.16im, 0.25 - 0.18im, 1.15 + 0.08im, -1.1 - 0.12im]

"""
    shifted_diagonal_exponential_nep([roots]; kwargs...)

Return a diagonal `SPMF_NEP` with entries `m_j(λ) = exp(λ - roots[j]) - 1`.
"""
function shifted_diagonal_exponential_nep(
    roots=SHIFTED_DIAGONAL_EXPONENTIAL_NEP_ROOTS;
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:clean,
    notes="Shifted diagonal exponential SPMF_NEP with roots away from the contour boundary."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    nep = SPMF_NEP(
        [
            Matrix(Diagonal(exp.(-root_values))),
            -Matrix{ComplexF64}(I, n, n),
        ],
        [S -> exp(S), S -> one(S)]
    )

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="shifted_diagonal_exponential_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const DIAGONAL_DELAY_LAMBERTW_NEP_ROOTS =
    ComplexF64[-0.22 + 0.14im, 0.18 - 0.2im, 1.15 + 0.08im, 0.95 - 0.72im]

"""
    diagonal_delay_lambertw_nep([roots]; kwargs...)

Return a diagonal `DEP` with chosen roots verified by construction.
"""
function diagonal_delay_lambertw_nep(
    roots=DIAGONAL_DELAY_LAMBERTW_NEP_ROOTS;
    τ=1.0,
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:clean,
    notes="Diagonal delay DEP with precomputed roots verified by -r + (r*exp(τ*r))*exp(-τ*r) = 0; equivalently r = W_0(τ*a)/τ on the principal branch."
)
    root_values = _complex_vector(roots)
    A = Matrix(Diagonal(root_values .* exp.(τ .* root_values)))
    nep = DEP([A], [τ])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="diagonal_delay_lambertw_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const SIMILARITY_DIAGONAL_NEP_ROOTS =
    ComplexF64[-0.31 + 0.11im, 0.26 - 0.19im, 1.22 + 0.1im, -1.18 - 0.16im]

const SIMILARITY_DIAGONAL_NEP_TRANSFORM = ComplexF64[
    1      0.2   -0.1   0.05
    0.1    1      0.15 -0.08
   -0.05   0.12   1      0.2
    0.08  -0.04   0.1    1
]

"""
    similarity_diagonal_nep([roots]; kwargs...)

Return a linear `PEP` similar to a diagonal known-spectrum problem.
"""
function similarity_diagonal_nep(
    roots=SIMILARITY_DIAGONAL_NEP_ROOTS;
    transform=SIMILARITY_DIAGONAL_NEP_TRANSFORM,
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:clean,
    notes="Similarity-transformed diagonal linear PEP with a fixed deterministic transform."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    Q = Matrix{ComplexF64}(transform)
    length(root_values) == size(Q, 1) == size(Q, 2) || error("transform size must match roots")
    A0 = -(Q * Matrix(Diagonal(root_values)) * inv(Q))
    nep = PEP([A0, Matrix{ComplexF64}(I, n, n)])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="similarity_diagonal_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const ILL_SCALED_DIAGONAL_NEP_ROOTS =
    ComplexF64[-0.33 + 0.13im, 0.21 - 0.24im, 1.18 + 0.12im, -1.16 - 0.18im]

const ILL_SCALED_DIAGONAL_NEP_SCALES = Float64[1e-6, 1e-2, 1e2, 1e6]

"""
    ill_scaled_diagonal_nep([roots]; kwargs...)

Return a diagonal linear `PEP` with controlled roots and varied component scales.
"""
function ill_scaled_diagonal_nep(
    roots=ILL_SCALED_DIAGONAL_NEP_ROOTS;
    scales=ILL_SCALED_DIAGONAL_NEP_SCALES,
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:diagnostic,
    notes="Ill-scaled diagonal linear PEP with roots controlled but diagonal component scales varying by orders of magnitude."
)
    root_values = _complex_vector(roots)
    scale_values = collect(scales)
    length(root_values) == length(scale_values) || error("scales length must match roots")
    nep = PEP([
        Matrix(Diagonal(-scale_values .* root_values)),
        Matrix(Diagonal(complex.(scale_values))),
    ])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="ill_scaled_diagonal_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const JORDAN_LIKE_NEP_ROOTS =
    ComplexF64[0.18 + 0.12im, 0.18 + 0.12im, -0.34 - 0.08im, 1.2 + 0.1im]

"""
    jordan_like_nep([roots]; kwargs...)

Return a linear `PEP` with a small Jordan block for a repeated eigenvalue.
"""
function jordan_like_nep(
    roots=JORDAN_LIKE_NEP_ROOTS;
    jordan_coupling=0.05 + 0.0im,
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:stress,
    notes="Jordan-like linear PEP with a repeated target root and a small off-diagonal coupling."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    J = Matrix(Diagonal(root_values))
    n >= 2 || error("jordan_like_nep requires at least two roots")
    J[1, 2] = jordan_coupling
    nep = PEP([-J, Matrix{ComplexF64}(I, n, n)])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="jordan_like_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end

const CLUSTERED_BOUNDARY_NEP_ROOTS =
    ComplexF64[0.62 + 0.02im, 0.625 - 0.015im, 0.67 + 0.01im, 0.675 - 0.02im]

"""
    clustered_boundary_nep([roots]; kwargs...)

Return a diagonal linear `PEP` with roots clustered near the contour boundary.
"""
function clustered_boundary_nep(
    roots=CLUSTERED_BOUNDARY_NEP_ROOTS;
    σ=0.0 + 0.0im,
    radius=0.65,
    solver_kwargs=(; N=256, k=3, K=2),
    expected_behavior=:stress,
    notes="Diagonal linear PEP with roots clustered close to the circular contour boundary."
)
    root_values = _complex_vector(roots)
    n = length(root_values)
    nep = PEP([Matrix(Diagonal(-root_values)), Matrix{ComplexF64}(I, n, n)])

    inside = [root for root in root_values if inside_ellipse(root; σ=σ, radius=radius)]
    outside = [root for root in root_values if !inside_ellipse(root; σ=σ, radius=radius)]

    return block_ss_known_spectrum_problem(;
        name="clustered_boundary_nep",
        nep=nep,
        inside=inside,
        outside=outside,
        σ=σ,
        radius=radius,
        solver_kwargs=solver_kwargs,
        expected_behavior=expected_behavior,
        notes=notes,
    )
end
