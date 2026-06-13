using LinearAlgebra
using Random

export SquareRegion, RectangularRegion, EllipseContour
export SIMResult, SIMContourDecision
export enclosing_contour, contour_parameters, inside_region, subdivide
export region_boundary, contour_boundary, collect_region_boundaries, collect_contour_boundaries
export sim_indicator, sim_screen, sim_screen_regions, sim_contour_decision

"""
    SquareRegion(center, half_side)

Axis-aligned square search region in the complex plane. The region is centered
at `center` and extends `half_side` in both the real and imaginary directions.
"""
struct SquareRegion{T<:Real}
    center::Complex{T}
    half_side::T
end

"""
    RectangularRegion(center, half_width, half_height)

Axis-aligned rectangular search region in the complex plane. The region is
centered at `center`, with half-width in the real direction and half-height in
the imaginary direction.
"""
struct RectangularRegion{T<:Real}
    center::Complex{T}
    half_width::T
    half_height::T
end

"""
    EllipseContour(center, radius)

Ellipse contour parameters compatible with contour integral solvers. A scalar
`radius` creates a circle; a length-two radius stores the real and imaginary
semi-axes.
"""
struct EllipseContour{T<:Real}
    center::Complex{T}
    radius::Tuple{T,T}
end

"""
    SIMResult

Result of SIM screening for one region. `indicator` is the Xi-Sun
projection-ratio screening value only; it is not an eigenvalue count.
"""
struct SIMResult{T<:Real,R}
    region::R
    indicator::T
    active::Bool
    threshold::T
    N::Int
    probe_norm::T
end

"""
    SIMContourDecision

Policy-layer recommendation from combining a `SIMResult` with optional contour
diagnostics. `recommended_kwargs` is reserved for future policy suggestions and
is currently empty.
"""
struct SIMContourDecision
    action::Symbol
    reason::Symbol
    recommended_kwargs::NamedTuple
end

_real_promote_type(center::Number, values::Real...) =
    promote_type(typeof(real(complex(center))), map(typeof, values)...)

function SquareRegion(center::Number, half_side::Real)
    half_side > 0 || error("half_side must be positive")
    T = _real_promote_type(center, half_side)
    SquareRegion{T}(Complex{T}(center), T(half_side))
end

function RectangularRegion(center::Number, half_width::Real, half_height::Real)
    half_width > 0 || error("half_width must be positive")
    half_height > 0 || error("half_height must be positive")
    T = _real_promote_type(center, half_width, half_height)
    RectangularRegion{T}(Complex{T}(center), T(half_width), T(half_height))
end

function EllipseContour(center::Number, radius::Real)
    radius > 0 || error("radius must be positive")
    T = _real_promote_type(center, radius)
    EllipseContour{T}(Complex{T}(center), (T(radius), T(radius)))
end

function EllipseContour(center::Number, radius)
    length(radius) == 2 || error("radius must be a scalar or a length-two collection")
    radius[1] > 0 || error("first radius must be positive")
    radius[2] > 0 || error("second radius must be positive")
    T = _real_promote_type(center, radius[1], radius[2])
    EllipseContour{T}(Complex{T}(center), (T(radius[1]), T(radius[2])))
end

_region_half_width(region::SquareRegion) = region.half_side
_region_half_height(region::SquareRegion) = region.half_side
_region_half_width(region::RectangularRegion) = region.half_width
_region_half_height(region::RectangularRegion) = region.half_height

"""
    enclosing_contour(region; shape=:circle)

Convert a square or rectangular region to an enclosing `EllipseContour`.
The default `:circle` is the circumscribed circle. `shape=:ellipse` returns an
axis-aligned ellipse with semi-axes equal to the region half-width and
half-height. SIM screening uses the circumscribed circle, so rectangular
screening can produce false positives that should later be filtered with
`inside_region`.
"""
function enclosing_contour(region::Union{SquareRegion,RectangularRegion}; shape=:circle)
    half_width = _region_half_width(region)
    half_height = _region_half_height(region)
    if shape == :circle
        return EllipseContour(region.center, sqrt(half_width^2 + half_height^2))
    elseif shape == :ellipse
        return EllipseContour(region.center, (half_width, half_height))
    else
        error("unknown enclosing contour shape: $shape")
    end
end

"""
    contour_parameters(contour)

Return named contour parameters `(σ=..., radius=...)` suitable for
`contour_beyn`, `contour_beyn_info`, `contour_block_SS`, and
`contour_block_SS_info`.
"""
contour_parameters(contour::EllipseContour) = (σ=contour.center, radius=contour.radius)

"""
    inside_region(region, z)

Return true if the complex point `z` lies inside or on the boundary of a square
or rectangular region.
"""
function inside_region(region::Union{SquareRegion,RectangularRegion}, z::Number)
    dz = complex(z) - region.center
    return abs(real(dz)) <= _region_half_width(region) &&
        abs(imag(dz)) <= _region_half_height(region)
end

"""
    subdivide(region)

Split a square or rectangular region into four equal axis-aligned child
regions.
"""
function subdivide(region::SquareRegion)
    child_half = region.half_side / 2
    offsets = (child_half + im*child_half,
               -child_half + im*child_half,
               child_half - im*child_half,
               -child_half - im*child_half)
    return [SquareRegion(region.center + offset, child_half) for offset in offsets]
end

function subdivide(region::RectangularRegion)
    child_half_width = region.half_width / 2
    child_half_height = region.half_height / 2
    offsets = (child_half_width + im*child_half_height,
               -child_half_width + im*child_half_height,
               child_half_width - im*child_half_height,
               -child_half_width - im*child_half_height)
    return [RectangularRegion(region.center + offset, child_half_width, child_half_height)
            for offset in offsets]
end

_boundary_tuple(z) = (x=real.(z), y=imag.(z), z=z)

"""
    region_boundary(region; n=2)

Return named coordinate vectors `(x=..., y=..., z=...)` for the closed
boundary of a square or rectangular region. `n` is the number of samples per
edge and must be at least 2.
"""
function region_boundary(region::Union{SquareRegion,RectangularRegion}; n::Integer=2)
    n >= 2 || error("n must be at least 2")
    half_width = _region_half_width(region)
    half_height = _region_half_height(region)
    x = collect(range(real(region.center) - half_width,
                      stop=real(region.center) + half_width,
                      length=n))
    y = collect(range(imag(region.center) - half_height,
                      stop=imag(region.center) + half_height,
                      length=n))
    z = [complex(xi, y[1]) for xi in x]
    append!(z, (complex(x[end], yi) for yi in y[2:end]))
    append!(z, (complex(xi, y[end]) for xi in reverse(x[1:end-1])))
    append!(z, (complex(x[1], yi) for yi in reverse(y[1:end-1])))
    return _boundary_tuple(z)
end

"""
    contour_boundary(contour; n=100, closed=true)
    contour_boundary(center, radius; n=100, closed=true)
    contour_boundary(; σ, radius, n=100, closed=true)

Return named coordinate vectors `(x=..., y=..., z=...)` for an ellipse or
circle contour. The two-argument and keyword forms match the contour-solver
convention of a center `σ` and scalar or length-two `radius`.
"""
function contour_boundary(contour::EllipseContour; n::Integer=100, closed::Bool=true)
    n > 0 || error("n must be positive")
    last_theta = closed ? 2*pi : 2*pi*(n - 1)/n
    theta = range(0, stop=last_theta, length=closed ? n + 1 : n)
    z = [contour.center + complex(contour.radius[1]*cos(t), contour.radius[2]*sin(t))
         for t in theta]
    return _boundary_tuple(z)
end

contour_boundary(center::Number, radius; params...) =
    contour_boundary(EllipseContour(center, radius); params...)

contour_boundary(; σ, radius, params...) =
    contour_boundary(σ, radius; params...)

_path_id(path::Tuple{}) = "root"
_path_id(path::Tuple) = join(path, ".")

function _validate_depth(depth)
    depth >= 0 || error("depth must be nonnegative")
    return Int(depth)
end

function _validate_path(path)
    all(i -> i isa Integer && 1 <= i <= 4, path) ||
        error("selected paths must contain child indices between 1 and 4")
    return Tuple(Int(i) for i in path)
end

function _normalize_selected_paths(selected)
    selected === nothing && return nothing
    if selected isa Tuple
        return [_validate_path(selected)]
    elseif selected isa AbstractVector && all(i -> i isa Integer, selected)
        return [_validate_path(selected)]
    elseif selected isa AbstractVector
        return [_validate_path(path) for path in selected]
    else
        error("selected must be a path tuple or a vector of path tuples")
    end
end

function _region_at_path(region, path::Tuple)
    current = region
    for child_index in path
        current = subdivide(current)[child_index]
    end
    return current
end

function _selected_prefixes(selected_paths)
    paths = Set{Tuple}()
    push!(paths, ())
    for path in selected_paths
        for i = 1:length(path)
            push!(paths, path[1:i])
        end
    end
    sorted_paths = collect(paths)
    sort!(sorted_paths, by=path -> (length(path), _path_id(path)))
    return sorted_paths
end

function _region_entries(region, depth::Integer, selected)
    selected_paths = _normalize_selected_paths(selected)
    if selected_paths !== nothing
        return [(region=_region_at_path(region, path), path=path)
                for path in _selected_prefixes(selected_paths)]
    end

    entries = []
    function visit(current_region, path)
        push!(entries, (region=current_region, path=path))
        length(path) == depth && return
        for (child_index, child) in enumerate(subdivide(current_region))
            visit(child, (path..., child_index))
        end
    end
    visit(region, ())
    return entries
end

function _region_metadata(region, path::Tuple)
    return (region=region,
            level=length(path),
            path=path,
            id=_path_id(path),
            parent_id=isempty(path) ? nothing : _path_id(path[1:end-1]),
            child_index=isempty(path) ? nothing : last(path))
end

"""
    collect_region_boundaries(root_region; depth=0, selected=nothing, n=2)

Return boundary coordinates and hierarchy metadata for a region subdivision.
If `selected` is supplied as a child-index path such as `(1, 4)`, or a vector
of paths such as `[(1,), (2, 3)]`, the root and all selected path prefixes are
returned.
"""
function collect_region_boundaries(
    root_region::Union{SquareRegion,RectangularRegion};
    depth::Integer=0,
    selected=nothing,
    n::Integer=2
)
    depth = _validate_depth(depth)
    return [begin
                boundary = region_boundary(entry.region; n=n)
                merge(_region_metadata(entry.region, entry.path), boundary)
            end for entry in _region_entries(root_region, depth, selected)]
end

"""
    collect_contour_boundaries(root_region; depth=0, selected=nothing, n=100, shape=:circle)

Return enclosing contour coordinates and hierarchy metadata for a region
subdivision. Contours are created with `enclosing_contour(region; shape=shape)`.
"""
function collect_contour_boundaries(
    root_region::Union{SquareRegion,RectangularRegion};
    depth::Integer=0,
    selected=nothing,
    n::Integer=100,
    shape=:circle
)
    depth = _validate_depth(depth)
    return [begin
                contour = enclosing_contour(entry.region; shape=shape)
                boundary = contour_boundary(contour; n=n)
                merge(merge(_region_metadata(entry.region, entry.path), (contour=contour,)), boundary)
            end for entry in _region_entries(root_region, depth, selected)]
end

function _unit_probe(::Type{T}, n::Integer; probe=nothing, seed=10) where {T<:Number}
    if probe === nothing
        rng = MersenneTwister(seed)
        v = randn(rng, real(T), n)
    else
        length(probe) == n || error("probe length must match the NEP size")
        v = Vector{T}(probe)
    end
    probe_norm = norm(v)
    probe_norm > 0 || error("probe vector must be nonzero")
    return Vector{T}(v / probe_norm), real(probe_norm)
end

function _sim_circle_sums(
    ::Type{T},
    nep::NEP,
    region::Union{SquareRegion,RectangularRegion};
    N::Integer,
    linsolvercreator,
    probe=nothing,
    seed=10
) where {T<:Number}
    N > 0 || error("N must be positive")
    iseven(N) || error("N must be even for the Xi-Sun half-grid SIM indicator")
    n = size(nep, 1)
    f, probe_norm = _unit_probe(T, n; probe=probe, seed=seed)
    contour = enclosing_contour(region; shape=:circle)
    center = contour.center
    radius = contour.radius[1]
    full_grid = zeros(T, n)
    half_grid = zeros(T, n)

    # Xi-Sun pmCIM indicator for a disk: I = norm(Pf_N ./ Pf_(N/2)) / sqrt(n),
    # where Pf_(N/2) uses the even-index half-grid {x_2, x_4, ..., x_N}.
    for j = 1:N
        theta = 2*pi*(j - 1)/N
        phase = exp(im*theta)
        z = center + radius*phase
        solver = create_linsolver(linsolvercreator, nep, z)
        term = radius * phase / N * lin_solve(solver, f)
        full_grid .+= term
        if iseven(j)
            half_grid .+= 2 .* term
        end
    end
    return full_grid, half_grid, probe_norm
end

function _projection_ratio_indicator(full_grid, half_grid)
    ratios = similar(full_grid)
    scale = max(norm(full_grid), norm(half_grid), 1)
    denom_tol = eps(real(eltype(full_grid))) * scale
    for i in eachindex(full_grid)
        if abs(half_grid[i]) <= denom_tol
            ratios[i] = abs(full_grid[i]) <= denom_tol ? zero(eltype(full_grid)) : complex(Inf)
        else
            ratios[i] = full_grid[i] / half_grid[i]
        end
    end
    return real(norm(ratios) / sqrt(length(full_grid)))
end

"""
    sim_indicator(nep, region; N=16, linsolvercreator=BackslashLinSolverCreator(), probe=nothing, seed=10)

Compute the Xi-Sun SIM projection-ratio indicator for active/inactive region
screening. The returned scalar is not an eigenvalue count and no eigenpairs are
computed. Use `sim_screen` to apply a threshold.
"""
function sim_indicator(
    nep::NEP,
    region::Union{SquareRegion,RectangularRegion};
    N::Integer=16,
    linsolvercreator=BackslashLinSolverCreator(),
    probe=nothing,
    seed=10
)
    full_grid, half_grid, _ = _sim_circle_sums(ComplexF64, nep, region;
                                               N=N,
                                               linsolvercreator=linsolvercreator,
                                               probe=probe,
                                               seed=seed)
    return _projection_ratio_indicator(full_grid, half_grid)
end

"""
    sim_screen(nep, region; N=16, threshold=0.1, ...)

Compute `sim_indicator` and return a `SIMResult` with `active =
indicator > threshold`.
"""
function sim_screen(
    nep::NEP,
    region::Union{SquareRegion,RectangularRegion};
    N::Integer=16,
    threshold::Real=0.1,
    linsolvercreator=BackslashLinSolverCreator(),
    probe=nothing,
    seed=10
)
    full_grid, half_grid, probe_norm = _sim_circle_sums(ComplexF64, nep, region;
                                                        N=N,
                                                        linsolvercreator=linsolvercreator,
                                                        probe=probe,
                                                        seed=seed)
    indicator = _projection_ratio_indicator(full_grid, half_grid)
    Tind = typeof(indicator)
    return SIMResult(region, indicator, indicator > threshold, Tind(threshold), Int(N), Tind(probe_norm))
end

"""
    sim_screen_regions(nep, regions; params...)

Apply `sim_screen` to each region and return a vector of `SIMResult`s.
"""
sim_screen_regions(nep::NEP, regions; params...) =
    [sim_screen(nep, region; params...) for region in regions]

function _info_field(info, name::Symbol, default)
    info === nothing && return default
    hasproperty(info, name) || return default
    return getproperty(info, name)
end

function _unclear_rank_gap(info)
    singular_values = _info_field(info, :singular_values, nothing)
    estimated_rank = _info_field(info, :estimated_rank, nothing)
    rank_drop_tol = _info_field(info, :rank_drop_tol, nothing)
    if singular_values === nothing || estimated_rank === nothing || rank_drop_tol === nothing
        return false
    end
    isempty(singular_values) && return false
    first_sv = singular_values[1]
    first_sv == 0 && return false
    r = Int(estimated_rank)
    1 <= r < length(singular_values) || return false
    normalized = singular_values ./ first_sv
    next_sv = normalized[r + 1]
    next_sv > rank_drop_tol / 10 && return true
    next_sv == 0 && return false
    return normalized[r] / next_sv < 10
end

"""
    sim_contour_decision(sim_result, contour_info=nothing; capacity_margin=1, residual_tol=nothing)

Policy helper that combines SIM screening with optional contour diagnostics.
SIM supplies only the active/inactive decision; rank, capacity, and residual
information are read only from contour diagnostics.
"""
function sim_contour_decision(sim_result::SIMResult, contour_info=nothing;
                              capacity_margin::Integer=1, residual_tol=nothing)
    if !sim_result.active
        return SIMContourDecision(:skip, :inactive, NamedTuple())
    end
    contour_info === nothing && return SIMContourDecision(:run_beyn, :active, NamedTuple())

    estimated_rank = _info_field(contour_info, :estimated_rank, 0)
    capacity = _info_field(contour_info, :capacity, 0)
    number_returned = _info_field(contour_info, :number_returned, nothing)
    residuals = _info_field(contour_info, :residuals, [])

    if number_returned !== nothing && number_returned == 0
        return SIMContourDecision(:subdivide_or_rerun, :active_no_eigenpairs, NamedTuple())
    end

    if residual_tol !== nothing && !isempty(residuals) && any(residuals .> residual_tol)
        return SIMContourDecision(:subdivide_or_rerun, :poor_residuals, NamedTuple())
    end

    if _unclear_rank_gap(contour_info)
        return SIMContourDecision(:subdivide_or_rerun, :unclear_rank_gap, NamedTuple())
    end

    if capacity > 0 && estimated_rank >= capacity - capacity_margin
        return SIMContourDecision(:run_block_SS_or_subdivide, :near_capacity_rank, NamedTuple())
    end

    return SIMContourDecision(:accept_beyn, :low_estimated_rank, NamedTuple())
end
