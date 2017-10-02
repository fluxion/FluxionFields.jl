import Base: convert

using FluxionArrays


# Abstract types


# Implements:
# `convert(FieldExpr, x)`
abstract type FieldExpr end

abstract type Dimension <: FieldExpr end
abstract type Field <: FieldExpr end
abstract type KnownDimension <: Dimension end


# Dimensions


struct UnknownDimension <: Dimension
    id :: Symbol

    UnknownDimension() = new(gensym("UD"))
end


struct BasicDimension <: KnownDimension
    id :: Symbol
    grid :: Grid

    BasicDimension(grid::Grid) = new(gensym("BD"), grid)
end


struct TransformedDimension <: KnownDimension
    id :: Symbol
    grid :: Grid
    base_dimension :: BasicDimension

    TransformedDimension(grid::Grid, base_dimension::BasicDimension) =
        new(gensym("TD"), grid, base_dimension)
end


# Stochastic metadata


# This object represents three possible states:
# a non-stochastic field expression, an unknown stochastic field expression,
# and a known stochastic field expression (which implies a known number of trajectories)
struct StochasticInfo
    trajectories :: Nullable{Int64}

    StochasticInfo() = new(nothing)
    StochasticInfo(trajectories::Int64) = new(trajectories)
end


struct StochasticMean <: FieldExpr
    fexpr :: FieldExpr
end


struct Sample <: FieldExpr
    fexpr :: FieldExpr
    trajectories :: Int64
end


struct StochasticDimension <: KnownDimension
end

const STOCHASTIC_DIMENSION = StochasticDimension()


# Fields


struct UnknownField <: Field
    id :: Symbol
    dimensions :: Array{Dimension, 1}
    stochastic :: Bool

    UnknownField(dimensions, stochastic=false) = new(gensym("UF"), dimensions, stochastic)
end


struct KnownField <: Field
    dimensions :: Array{KnownDimension, 1}
    data :: BackendArray
    stochastic_info :: StochasticInfo

    KnownField(dimensions, data) = new(dimensions, data, StochasticInfo())
    KnownField(dimensions, data, stochastic_info) = new(dimensions, data, stochastic_info)
end


struct FieldExprNode <: FieldExpr
    head :: Symbol
    args :: Array{FieldExpr, 1}
end


struct Fourier <: FieldExpr
    fexpr :: FieldExpr
    inverse :: Bool
    dimensions :: Dict{Dimension, Dimension}
end


struct Diff <: FieldExpr
    field :: UnknownField
    dimensions :: Dict{Dimension, Int}
end


struct RandomVariable <: FieldExpr
    distribution :: Symbol
    dimensions :: Array{Dimension, 1}
end


# Methods


grid(dim::BasicDimension) = dim.grid
grid(dim::TransformedDimension) = dim.grid


convert(::Type{FieldExpr}, val::Number) = KnownField([], to_backend(reshape([val], ())))
convert(::Type{FieldExpr}, val::FieldExpr) = val

to_subnodes(fexpr::FieldExpr) = nothing
to_subnodes(fexpr::FieldExprNode) = fexpr.args
to_subnodes(fexpr::Fourier) = [fexpr.fexpr]
to_subnodes(fexpr::Diff) = [fexpr.field]
to_subnodes(fexpr::StochasticMean) = [fexpr.fexpr]

from_subnodes(fexpr::FieldExprNode, subnodes) = FieldExprNode(fexpr.head, subnodes)
from_subnodes(fexpr::Fourier, subnodes) = Fourier(subnodes[1], fexpr.inverse, fexpr.dimensions)
from_subnodes(fexpr::Diff, subnodes) = Diff(subnodes[1], fexpr.dimensions)
from_subnodes(fexpr::StochasticMean, subnodes) = StochasticMean(subnodes[1])

known_dimensions(x) = []
known_dimensions(x::UnknownField) = [dim for dim in x.dimensions if isa(dim, KnownDimension)]
known_dimensions(x::KnownField) = x.dimensions
known_dimensions(x::KnownDimension) = [x]
known_dimensions(x::TransformedDimension) = [x]

unknown_dimensions(x::KnownField) = []
unknown_dimensions(x::UnknownField) = [dim for dim in x.dimensions if isa(dim, UnknownDimension)]
unknown_dimensions(x::KnownDimension) = []
unknown_dimensions(x::TransformedDimension) = []
unknown_dimensions(x::UnknownDimension) = [x]


function known_dimensions_expr(fexpr)
    inspect(PersistentSet{Dimension}(), fexpr) do controls, dims, x
        if typeof(fexpr) == Fourier
            skip_children!(controls)
            inner_dims = known_dimensions_expr(fexpr.fexpr)
            union(dims, [get(fexpr.dimensions, dim, dim) for dim in inner_dims])
        else
            union(dims, known_dimensions(x))
        end
    end
end


function unknown_dimensions_expr(fexpr)
    inspect(PersistentSet{Dimension}(), fexpr) do controls, dims, x
        union(dims, unknown_dimensions(x))
    end
end


function remove_unknown_dimensions(f::UnknownField)
    remove_dimensions(f, unknown_dimensions(f)...)
end


function remove_dimensions(f::UnknownField, dims...)
   to_remove = Set{Dimension}(dims)
   UnknownField([dim for dim in f.dimensions if !in(dim, to_remove)])
end


function replace_expr(fexpr::FieldExpr, replacements...)
    modify(fexpr, Dict(replacements...)) do controls, x, replacement_dict
        get(replacement_dict, x, x)
    end
end
