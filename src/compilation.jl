using FunctionalCollections
using FluxionArrays


function rearrange_dims(arr, arr_dims, target_dims)
    """
    Contracts:
    - elements of `arr_dims` and `target_dims` support the `==` operator
    - all elements of `arr_dims` are different
    - each element of `arr_dims` is present in `target_dims` exactly once
    - `arr` needs to support `reshape()`, `permutedims()` and `size()` methods.
    """
    if length(target_dims) == 0 || length(arr_dims) == 0
        return arr
    end

    old_shape = [size(arr, dim_num) for dim_num in 1:length(arr_dims)]

    # transpose the array so that the order of its dimensions
    # is the same as in `target_dims`
    dims_correspondence = [
        (findfirst(target_dims .== dim), old_dim_num)
        for (old_dim_num, dim) in enumerate(arr_dims)]
    new_dims_order = [old_dim_num for (new_dim_num, old_dim_num) in sort(dims_correspondence)]

    # Don't permute if the required order is already achieved
    if new_dims_order != collect(1:length(new_dims_order))
        arr = permutedims(arr, new_dims_order)
        new_shape = old_shape[new_dims_order]
        new_dims = arr_dims[new_dims_order]
    else
        new_shape = old_shape
        new_dims = arr_dims
    end

    # We may need to reshape the array adding new array dimensions (of size 1)
    # in place of unused array dimensions.
    # Currently just checking for the difference in shape lengths;
    # we probably don't need to add dummy dimensions if broadcasting works for us,
    # but that's harder to check
    if length(new_shape) != length(target_dims)

        new_dim_num = 1
        target_shape = Array{Any, 1}(length(target_dims))
        for target_dim_num = 1:length(target_dims)
            if (new_dim_num <= length(new_dims)
                    && target_dims[target_dim_num] == new_dims[new_dim_num])
                target_shape[target_dim_num] = new_shape[new_dim_num]
                new_dim_num += 1
            else
                target_shape[target_dim_num] = 1
            end
        end

        arr = reshape(arr, target_shape...)
    end

    return arr
end


function ensure_shape(arr, arr_shape, target_shape)

    if arr_shape == target_shape
        return arr
    end

    if length(arr_shape) < length(target_shape)
        extend_dims = length(target_shape) - length(arr_shape)
        extend_shape = [1 for i in 1:extend_dims]
        arr_shape = [arr_shape..., extend_shape...]
        arr = reshape(arr, arr_shape...)
    end

    outer = [
        (arr_dim == 1 && target_dim > 1) ? target_dim : 1
        for (arr_dim, target_dim) in zip(arr_shape, target_shape)
    ]

    repeat(arr, outer=outer)
end


function map_accum{T}(func, state, vals::Array{T, 1}, args...)
    new_vals = Array{Any, 1}(length(vals))
    for (i, val) in enumerate(vals)
        state, val = func(state, val, args...)
        new_vals[i] = val
    end
    state, new_vals
end


immutable ArrayExprBuildingState
    unknowns :: PersistentHashMap{FieldExpr, Unknown}
    random_states :: Nullable{Pair{Unknown, ArrayExpr}}
    random_variables :: PersistentHashMap{RandomVariable, ArrayExpr}
end

ArrayExprBuildingState() = ArrayExprBuildingState(
    PersistentHashMap{FieldExpr, Unknown}(),
    nothing,
    PersistentHashMap{RandomVariable, ArrayExpr}())


function add_if_missing(gen_val, state::ArrayExprBuildingState, key::FieldExpr)
    if haskey(state.unknowns, key)
        state, state.unknowns[key]
    else
        val = gen_val(key)
        new_unknowns = assoc(state.unknowns, key, val)
        new_state = ArrayExprBuildingState(
            new_unknowns, state.random_states, state.random_variables)
        new_state, val
    end
end


function set_random_state(state::ArrayExprBuildingState, rs::Unknown)
    ArrayExprBuildingState(
        state.unknowns, Pair{Unknown, ArrayExpr}(rs, rs), state.random_variables)
end


function update_random_state(state::ArrayExprBuildingState, rs::ArrayExpr)
    input_rs, output_rs = get(state.random_states)
    ArrayExprBuildingState(
        state.unknowns, Pair{Unknown, ArrayExpr}(input_rs, rs), state.random_variables)
end


function add_random_variable(state::ArrayExprBuildingState, fexpr::RandomVariable, aexpr::ArrayExpr)
    ArrayExprBuildingState(
        state.unknowns, state.random_states,
        assoc(state.random_variables, fexpr, aexpr))
end

function clear_random_variables(state::ArrayExprBuildingState)
    ArrayExprBuildingState(
        state.unknowns, state.random_states,
        PersistentHashMap{RandomVariable, ArrayExpr}())
end

function set_random_variables(state::ArrayExprBuildingState, random_variables)
    ArrayExprBuildingState(
        state.unknowns, state.random_states,
        random_variables)
end



random_states_set(state::ArrayExprBuildingState) = !isnull(state.random_states)


random_states(state::ArrayExprBuildingState) = get(state.random_states)
output_random_state(state::ArrayExprBuildingState) = random_states(state)[2]


function unknown_array_like(x::UnknownDimension)
    return unknown_array(shape=[], eltype=Float64)
end

function unknown_array_like(x::UnknownField)
    return unknown_array(
        shape=[Nullable{Int64}(length(vertices(grid(dim)))) for dim in known_dimensions(x)],
        eltype=Float64)
end


_as_array_expr(state, x::Number, target_dimensions) = state, lazy_array(x)
_as_array_expr(state, x::KnownDimension, target_dimensions) =
    state, rearrange_dims(lazy_array(vertices(grid(x))), [x], target_dimensions)
_as_array_expr(state, x::TransformedDimension, target_dimensions) =
    state, rearrange_dims(lazy_array(vertices(grid(x))), [x], target_dimensions)
_as_array_expr(state, x::KnownField, target_dimensions) =
    state, rearrange_dims(x.data, known_dimensions(x), target_dimensions)
_as_array_expr(state, x::UnknownDimension, target_dimensions) =
    add_if_missing(unknown_array_like, state, x)

function _as_array_expr(state, x::UnknownField, target_dimensions)
    state, arr_var = add_if_missing(unknown_array_like, state, x)
    state, rearrange_dims(arr_var, known_dimensions(x), target_dimensions)
end

function _as_array_expr(state, x::Fourier, target_dimensions)

    dims_dict = map(reverse, x.dimensions)

    fft_axes = []
    target_dimensions_inner = []
    for (i, dim) in enumerate(target_dimensions)
        if haskey(dims_dict, dim)
            push!(fft_axes, i)
            push!(target_dimensions_inner, dims_dict[dim])
        else
            push!(target_dimensions_inner, dim)
        end
    end

    state, expr = _as_array_expr(state, x.fexpr, target_dimensions_inner)

    if x.inverse
        result = ifft(expr, fft_axes)
    else
        result = fft(expr, fft_axes)
    end
    state, result
end


function _as_array_expr(state, fexpr::Sample, target_dimensions)
    target_dimensions = [
        x === STOCHASTIC_DIMENSION ? KnownStochasticDimension(fexpr.trajectories) : x
        for x in target_dimensions]

    # Random variables inside a sample() call are independent from the rest
    outer_rvs = state.random_variables
    inner_state = clear_random_variables(state)
    state, aexpr = _as_array_expr(inner_state, fexpr.fexpr, target_dimensions)
    state = set_random_variables(state, outer_rvs)
    state, aexpr
end


function _as_array_expr(state, fexpr::StochasticMean, target_dimensions)
    trs = trajectories_number(fexpr.fexpr)
    target_dimensions = [target_dimensions..., KnownStochasticDimension(trs)]
    state, aexpr = _as_array_expr(state, fexpr.fexpr, target_dimensions)
    state, sum(aexpr, length(target_dimensions)) ./ trs
end


immutable KnownStochasticDimension <: KnownDimension
    id :: Symbol
    trajectories :: Int64
    KnownStochasticDimension(trajectories) = new(gensym(), trajectories)
end


function _as_array_expr(state, fexpr::RandomVariable, target_dimensions)

    if haskey(state.random_variables, fexpr)
        return state, state.random_variables[fexpr]
    end

    if fexpr.distribution == :random_normal
        sdims = find(dim -> typeof(dim) == KnownStochasticDimension, target_dimensions)
        @assert length(sdims) == 1
        random_dimensions = [fexpr.dimensions..., target_dimensions[sdims[1]]]

        shape = [
            typeof(x) == KnownStochasticDimension ? x.trajectories : length(vertices(grid(x)))
            for x in random_dimensions]

        if random_states_set(state)
            rs = output_random_state(state)
        else
            rs = unknown_random_state()
            state = set_random_state(state, rs)
        end

        # TODO: we need to either use the known mean/std here,
        # or expose the basic random_normal(rs, shape...) from Arrays.
        res, new_rs = FluxionArrays.random_normal(rs, 0, 1, shape...)

        res = rearrange_dims(res, random_dimensions, target_dimensions)

        state = update_random_state(state, new_rs)
        state = add_random_variable(state, fexpr, res)

        return state, res
    end

end


function _as_array_expr(state, fexpr::FieldExprNode, target_dimensions)
    const broadcasted_operators = Set([:+, :-, :*, :/, :^])

    op = fexpr.head
    state, args = map_accum(_as_array_expr, state, fexpr.args, target_dimensions)
    if in(op, broadcasted_operators)
        res = broadcast(getfield(Base, op), args...)
    else
        res = getfield(Base, op)(args...)
    end
    state, res
end


function as_array_expr(state, fexpr::FieldExpr, target_dimensions)

    state, aexpr = _as_array_expr(state, fexpr, target_dimensions)

    if is_stochastic(fexpr)
        tr_num = trajectories_number(fexpr)
    else
        tr_num = nothing
    end
    target_shape = [
        dim === STOCHASTIC_DIMENSION ? tr_num : length(vertices(grid(dim)))
        for dim in target_dimensions]

    shape = size(aexpr)
    @assert !isa(shape, ArrayExpr)
    shape = collect(shape)

    aexpr = ensure_shape(aexpr, shape, target_shape)

    state, aexpr
end


immutable WrappedArrayFunction
    array_function
    uses_random_state :: Bool
end


function wrapped_array_function(returns, parameters)

    array_returns = []

    state = ArrayExprBuildingState()

    for (fexpr, result_dimensions) in returns
        state, aexpr = as_array_expr(state, fexpr, result_dimensions)
        push!(array_returns, aexpr)
    end

    array_parameters = [
        get(state.unknowns, parameter, unknown_array())
        for parameter in parameters]

    uses_random_state = !isnull(state.random_states)

    if uses_random_state
        input_rs, output_rs = get(state.random_states)
        array_parameters = [input_rs, array_parameters...]
        array_returns = [output_rs, array_returns...]
    end

    afunc = array_function(array_returns, array_parameters)

    WrappedArrayFunction(afunc, uses_random_state)
end


immutable WrappedBackendArrayFunction
    backend_array_function
    uses_random_state :: Bool
end


function compile_wrapped_array_function(backend::ArrayBackend, afunc::WrappedArrayFunction)
    WrappedBackendArrayFunction(
        compile_array_function(backend, afunc.array_function), afunc.uses_random_state)
end


import FluxionArrays: array_backend

array_backend(wbafunc::WrappedBackendArrayFunction) = array_backend(wbafunc.backend_array_function)


function (wcafunc::WrappedBackendArrayFunction)(args...; seed=nothing)
    if wcafunc.uses_random_state
        backend = array_backend(wcafunc.backend_array_function)
        rs = random_state(backend, seed=seed)
        rs, result = wcafunc.backend_array_function(rs, args...)
        return result
    else
        wcafunc.backend_array_function(args...)
    end
end


immutable FieldFunction
    wrapped_array_function
    output_dimensions
end


function field_function(returns, parameters)
    wafunc = wrapped_array_function(returns, parameters)
    output_dimensions = [rd for (_, rd) in returns]
    FieldFunction(wafunc, output_dimensions)
end


immutable BackendFieldFunction
    wrapped_backend_array_function :: WrappedBackendArrayFunction
    output_dimensions
end


function compile_field_function(backend::ArrayBackend, ffunc::FieldFunction)
    BackendFieldFunction(
        compile_wrapped_array_function(backend, ffunc.wrapped_array_function),
        ffunc.output_dimensions)
end


unwrap_field(x::Number) = x
unwrap_field(x::KnownField) = x.data

function (bffunc::BackendFieldFunction)(args...; seed=nothing)

    args = map(unwrap_field, args)

    if bffunc.wrapped_backend_array_function.uses_random_state
        backend = array_backend(bffunc.wrapped_backend_array_function)
        rs = random_state(backend, seed=seed)
        results = bcfunc.wrapped_backend_array_function(rs, args...)
        results = results[2:end]
    else
        results = bffunc.wrapped_backend_array_function(args...)
        if length(bffunc.output_dimensions) == 1
            results = [results]
        end
    end

    wrapped_results = []
    for (result, odims) in zip(results, bffunc.output_dimensions)
        push!(wrapped_results, known_field(result, odims))
    end

    if length(wrapped_results) == 1
        return wrapped_results[1]
    else
        return tuple(wrapped_results)
    end
end
