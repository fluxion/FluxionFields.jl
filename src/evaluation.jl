function evaluate_to_array(backend::ArrayBackend, x::Number, target_dimensions; seed=nothing)
    evaluate(backend, convert(FieldExpr, x), target_dimensions)
end

function evaluate_to_array(backend::ArrayBackend, fexpr::FieldExpr, target_dimensions; seed=nothing)
    af = wrapped_array_function([fexpr => target_dimensions], [])
    cf = compile_wrapped_array_function(backend, af)
    cf(seed=seed)
end


function evaluate(backend::ArrayBackend, fexpr::FieldExpr, target_dimensions=nothing; seed=nothing)

    if target_dimensions === nothing
        target_dimensions = collect(known_dimensions_expr(fexpr))
    end

    if !is_stochastic(fexpr)
        dims = target_dimensions
    else
        dims = [target_dimensions..., STOCHASTIC_DIMENSION]
    end

    ffunc = field_function([fexpr => dims], [])
    cffunc = compile_field_function(backend, ffunc)
    cffunc()
end

evaluate(args...; kwds...) = evaluate(default_array_backend(), args...; kwds...)



function known_field{T <: KnownDimension}(x::BackendArray, dimensions::Array{T, 1})
    @assert ndims(x) == length(dimensions)
    sdims = find(dimensions .== STOCHASTIC_DIMENSION)
    if length(sdims) == 0
        regular_dims = dimensions
        stochastic_info = StochasticInfo()
    else
        regular_dims = [dimensions[1:sdims[1]-1]..., dimensions[sdims[1]+1:end]...]
        stochastic_info = StochasticInfo(size(x, sdims[1]))
    end
    KnownField(regular_dims, x, stochastic_info)
end
