function ordered_dimensions(fexpr::FieldExpr, preferred_dim_order::Array{Dimension, 1})

    used_dims = known_dimensions_expr(fexpr)
    preferred_dims = Set(preferred_dim_order)

    used_dims_not_in_pd = setdiff(used_dims, preferred_dims)

    [
        (dim for dim in preferred_dim_order if dim in used_dims)...,
        (dim for dim in used_dims_not_in_pd)...]
end


function join_fields(fs::Array{KnownField, 1}, dim::KnownDimension)

    @assert length(fs) == length(vertices(grid(dim)))

    # assuming that all fields in `fs` have the same dimensions!
    result_dims = [known_dimensions(fs[1])..., dim]
    result_dim_sizes = [length(vertices(grid(dim))) for dim in result_dims]

    backend = array_backend(fs[1].data) # FIXME
    arrs = [evaluate_to_array(backend, f, known_dimensions(fs[1])) for f in fs]
    arrs_j = [from_backend(arr) for arr in arrs]

    result_data_j = Array{eltype(arrs[1])}(result_dim_sizes...)
    colons = repeat([:], outer=length(result_dim_sizes)-1)
    for (i, arr) in enumerate(arrs_j)
        setindex!(result_data_j, arr, colons..., i)
    end
    result_data = to_backend(backend, result_data_j)

    known_field(result_data, result_dims)
end


function join_fields(
        fields::Array{KnownField, 1},
        new_dimension::KnownDimension, new_dimension_grid, generic_field)

    # FIXME: setting uniform=True for the time being to make plotting easier
    new_dim = TransverseDimension(new_dimension.name, new_dimension_grid, uniform=True)

    # Currently we're just attaching the new dimension in front,
    # but it is possible to preserve the position it has in the original field.
    new_field = KnownField([generic_field.dimensions..., new_dim])

    for (i, field) in enumerate(fields)
        # FIXME: this will require support in the array layer
        new_field.data[i] = field.data
    end

    return new_field
end


# Currently replacing with FFTs
# can be also replaced using finite differences.

function replace_differentials(fexpr::FieldExprNode)
    FieldExprNode(fexpr.head, [replace_differentials(arg) for arg in fexpr.args])
end

function replace_differentials(fexpr::Diff)
    field = fexpr.field

    result = convert(FieldExpr, field)

    dims = collect(keys(fexpr.dimensions))
    fdims = [fourier_space(dim) for dim in dims]

    forward_pairs = Dict(dim => fdim for (dim, fdim) in zip(dims, fdims))
    inverse_pairs = Dict(fdim => dim for (dim, fdim) in zip(dims, fdims))

    result = fourier(result, forward_pairs...)
    result = convert(FieldExpr, result)

    for (dim, order) in fexpr.dimensions
        result = result * (2pi * 1im * forward_pairs[dim])^order
    end
    result = fourier(result, inverse_pairs..., inverse=true)
    result = convert(FieldExpr, result)

    return result
end

function replace_differentials(fexpr)
    fexpr
end


fftfreq(n::Integer, d::Real) = [(0:(n-1)>>1)..., (-(n>>1):-1)...] ./ (d * n)

function fourier_space(x::KnownDimension)
    # possible parameters: normalization, ordered/unordered frequencies, discrete/continuous
    # (that `integral(f(x), x)` == either `integral(f(k), k)` or `sum(f(k), k)`)

    gx = grid(x)

    @assert typeof(gx) == RegularGrid

    freqs = fftfreq(gx.points, (gx.right - gx.left) / gx.points)
    TransformedDimension(ArbitraryGrid(freqs), x)
end


function has_free_random_variables(fexpr::FieldExpr)
    inspect(false, fexpr) do controls, free_rvs, x
        if typeof(x) == Sample
            skip_children!(controls)
        elseif !free_rvs && typeof(x) == RandomVariable
            free_rvs = true
        end

        free_rvs
    end
end


function _trajectories_number(fexpr::FieldExpr)
    # Perhaps can be replaced by a special field in the nodes and
    # the propagation of the stochastic info.
    inspect(Nullable{Int64}(), fexpr) do controls, trs, x
        if typeof(x) == StochasticMean
            skip_children!(controls)
        elseif typeof(x) == KnownField && !isnull(x.stochastic_info.trajectories)
            trs = x.stochastic_info.trajectories
        elseif typeof(x) == Sample
            trs = Nullable{Int64}(x.trajectories)
        end

        trs
    end
end


trajectories_number(fexpr) = get(_trajectories_number(fexpr))
is_stochastic(fexpr) = !isnull(_trajectories_number(fexpr))
