__precompile__()

module FluxionFields

include("grids.jl")
export RegularGrid, ArbitraryGrid
export vertices, is_monotonic, grid_step

include("types.jl")
export is_leaf, leaf_val, grid
export KnownField, KnownDimension, BasicDimension, UnknownDimension, UnknownField, Dimension, FieldExpr
export BasicDimension
export remove_unknown_dimensions, remove_dimensions
export known_dimensions, unknown_dimensions, replace_expr
export fourier_space
export STOCHASTIC_DIMENSION

include("compilation.jl")
export field_function
export compile_wrapped_array_function
export wrapped_array_function
export compile_field_function

include("helpers.jl")
export join_fields, replace_differentials
export find_dimension_order
export has_free_random_variables, is_stochastic, trajectories_number

include("evaluation.jl")
export evaluate
export evaluate_to_array
export known_field

include("operators.jl")
export Diff, diff, Equation, equal, fourier, sample
export random_normal
export stochastic_mean

include("show.jl")

include("traversal.jl")
export inspect, modify, inspect_and_modify, skip_children!

end # module
