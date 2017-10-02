import Base: diff


immutable Equation
    left :: FieldExpr
    right :: FieldExpr
end


for op in (:+, :-, :*, :/, :^, :equal)
    if op == :equal
        @eval ($op)(x::FieldExpr, y::FieldExpr) = Equation(x, y)
    else
        @eval begin
            @eval import Base: $op
            ($op)(x::FieldExpr, y::FieldExpr) = FieldExprNode($(Meta.quot(op)), [x, y])
        end
    end

    # ^ operator is a special case since there is a shortcut for ^(x, y::Integer) in Base,
    # that interferes with our generic definitions.
    # So we have to add specific definitions for Integer to avoid ambiguity.

    for ext_type in (Number, (op == :^ ? (Integer,) : ())...)
        @eval begin
            ($op)(x::FieldExpr, y::$ext_type) = ($op)(x, convert(FieldExpr, y))
            ($op)(x::$ext_type, y::FieldExpr) = ($op)(convert(FieldExpr, x), y)
        end
    end
end


for builtin_function in (:abs, :cosh)
    @eval begin
        import Base: $builtin_function

        function $builtin_function(x::FieldExpr, args...)
            # This check is not necessary, [x, args...] should work anyway,
            # but it currently doesn't because of bug #17003
            if length(args) > 0
                args = [convert(FieldExpr, arg) for arg in args]
                args = [x, args...]
            else
                args = [x]
            end
            FieldExprNode($(Meta.quot(builtin_function)), args)
        end
    end
end


# Taken from StatsBase.jl
function countmap{T}(x::AbstractArray{T,1})
    cm = Dict{T,Int}()
    for v in x
        cm[v] = get(cm, v, 0) + 1
    end
    return cm
end


function diff(f::UnknownField, dimensions...)
    # here we can check that the field actually have the required dimensions
    grouped_dims = countmap(collect(dimensions))
    Diff(f, grouped_dims)
end


function fourier(fexpr::FieldExpr, dimension::Dimension; inverse=false)
    fourier(fexpr, Set{Dimension}([dimension]), inverse=inverse)
end


function fourier(fexpr::FieldExpr, dimension_pairs...; inverse=false)
    if !inverse
        for (dim, fdim) in dimension_pairs
            @assert fdim.base_dimension === dim
        end
        dimension_dict = Dict{Dimension, Dimension}(dimension_pairs)
        Fourier(fexpr, false, dimension_dict)
    else
        for (fdim, dim) in dimension_pairs
            @assert fdim.base_dimension === dim
        end
        dimension_dict = Dict{Dimension, Dimension}(dimension_pairs)
        Fourier(fexpr, true, dimension_dict)
    end
end


function stochastic_mean(fexpr::FieldExpr)
# If the field is not stochastic, should we throw an error?
    if is_stochastic(fexpr)
        StochasticMean(fexpr)
    else
        error("Requires a stochastic field")
    end
end


function stochastic_std(fexpr::FieldExpr)
# If the field is not stochastic, should we throw an error?
    if is_stochastic(fexpr)
        sqrt(stochastic_mean(fexpr^2) - stochastic_mean(fexpr)^2)
    else
        error("Requires a stochastic field")
    end
end


function sample(fexpr::Number, trajectories)
    return Sample(convert(FieldExpr, fexpr), trajectories)
end

function sample(fexpr::FieldExpr, trajectories)
    return Sample(fexpr, trajectories)
end


function random_normal(dimensions, mean, std)
    # TODO: here we need to add the dV factor too, depending on what dimensions are requested
    RandomVariable(:random_normal, dimensions) * std + mean
end

