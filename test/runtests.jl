using FluxionArrays
using FluxionFields
using Base.Test


eval_get(args...; kwds...) = from_backend(
    evaluate_to_array(default_array_backend(), args...; kwds...))


@testset "Fields" begin


@testset "stage1" begin
    # That's what happens at the start of the integration when the initial condition
    # has to be evaluated into an array.

    x = BasicDimension(RegularGrid(0, 1, 20))
    y = BasicDimension(RegularGrid(0, 1, 30))

    x_arr = vertices(grid(x))
    y_arr = vertices(grid(y))

    f = x + y + 1

    # Then there would be some call like `integrate(eqn, f, ...)`
    # and inside it we will have to obtain the initial array representing `f`.

    f_arr = from_backend(evaluate(f, [x, y]).data)
    @test f_arr == x_arr .+ y_arr' .+ 1

    f_arr = eval_get(f, [x, y])
    @test f_arr == x_arr .+ y_arr' .+ 1

end


@testset "stage2" begin

    # That's what happens inside the integration.
    # An array function is prepared and then evaluated on each step.

    x = BasicDimension(RegularGrid(0, 1, 20))
    y = BasicDimension(RegularGrid(0, 1, 30))
    t = UnknownDimension()
    f = UnknownField([t, x, y])

    f_initial = x + y + 1
    right_side = f + t + x + 2

    # Variant 1: we treat known and unknown dimensions differently,
    # so f(t, x, y) is understood as an unknown array of two (known) dimensions, x and y.
    #
    # right_side_f = as_array_function(right_side, [f, t], [x, y])
    #
    # We need to ensure that the return value of this call can be supplied to the next call.
    # Currently it seems that the order of known dimensions in `f` is fixed
    # (since it is a `Field`, not `FieldExpr`), so all the information is available.
    #
    # What if we want to make a function that takes an array of `t` and a corresponding known
    # field `f(t, x, y)`? How do we distinguish between these cases? Do we actually need this case?

    # Variant 2: known and unknown dimensions are the same. E.g. an unknown field
    # for finding a stationary state of a heat equation has no unknown dimensions.
    # So every solver will do its own postprocessing of the equation,
    # replacing the field if necessary.
    f_step = remove_unknown_dimensions(f)
    @test known_dimensions(f_step) == [x, y]
    @test unknown_dimensions(f_step) == []
    right_side_step = replace_expr(right_side, f => f_step)
    # Do we automatically assume that an unknown dimension is a scalar?
    right_side_f = compile_wrapped_array_function(
        default_array_backend(),
        wrapped_array_function([right_side_step => [x, y]], [f_step, t]))

    # Then we can call the function as
    f_arr = eval_get(f_initial, [x, y])
    t = 3
    new_f_arr = from_backend(right_side_f(f_arr, t))

    x_arr = vertices(grid(x))
    y_arr = vertices(grid(y))

    @test new_f_arr == x_arr .+ y_arr' .+ 1 .+ t .+ x_arr .+ 2

end


@testset "stage3" begin

    # Finally, this is what's happening inside a sampler
    # We have an array that we need to wrap in a field
    # and present to the user for him to do operations on

    x = BasicDimension(RegularGrid(0, 1, 20))
    y = BasicDimension(RegularGrid(0, 1, 30))
    t = BasicDimension(RegularGrid(4, 5, 10))

    x_arr = vertices(grid(x))
    y_arr = vertices(grid(y))
    t_arr = vertices(grid(t))

    f_arr = eval_get(x + y + 1, [x, y])

    f = known_field(to_backend(f_arr), [x, y]) # this is passed to the user

    results = Array{KnownField, 1}()

    # defined by the user
    sampler(f_val, t_val) = f_val + x + t_val

    # constructed by the integrator
    f_unknown = UnknownField([x, y])
    t_unknown = UnknownDimension()
    sampler_expr = sampler(f_unknown, t_unknown)

    # The function will be a part of Fluxion, it's too high-level for this package
    # observable_dims = ordered_dimensions(sampler_expr, [x, y])
    observable_dims = [x, y]

    sampler_func = field_function([sampler_expr => observable_dims], [f_unknown, t_unknown])
    sampler_func = compile_field_function(default_array_backend(), sampler_func)

    # The integration loop
    for t_val in vertices(grid(t))
        observable = sampler_func(f, t_val)
        push!(results, observable)
    end

    samples = join_fields(results, t)

    @test known_dimensions(samples) == [x, y, t]

    samples_data = eval_get(samples, [x, y, t])
    reference_data = (
        reshape(x_arr, (length(x_arr), 1, 1)) .* 2
        .+ reshape(y_arr, (1, length(y_arr), 1))
        .+ 1
        .+ reshape(t_arr, (1, 1, length(t_arr))))

    @test isapprox(samples_data, reference_data)

end


#=
@testset "regular randoms" begin

    # Typical example: adding Wigner noise
    x = BasicDimension(RegularGrid(-5, 5, 64))
    classic_state = 1 / cosh(10 - x)

    k = fourier_space(x)
    classic_state_k = fourier(classic_state, x => k)
    wigner_state_k = classic_state_k + FluxionFields.random_normal([k], 0, 1/2) # TODO: should be complex
    wigner_state = fourier(wigner_state_k, k => x, inverse=true)

    # in the integrate()
    # we have the number of trajectories and a seed that we need to use to instantiate the array

    # Trajectories are passed to sample(), because there may be an expression which has
    # different number of trajectories in subtrees, e.g.
    # mean(sample(expr1, trajectories=100)) + mean(sample(expr2, trajectories=200))
    # But they also can be passed to as_array()

    # Seed is passed to as_array, because we will have to create a random state internally
    # to sample the randoms, and creating random states inside the generated function
    # will be technically difficult.

    wigner_state_sampled = sample(wigner_state, 100)
    arr1 = as_array(wigner_state_sampled, [x, STOCHASTIC_DIMENSION], seed=123)
    @test size(arr1) == (64, 100)

    arr2 = as_array(wigner_state_sampled, [x, STOCHASTIC_DIMENSION], seed=123) # same array
    @test isapprox(arr1, arr2)

    arr3 = as_array(wigner_state_sampled, [x, STOCHASTIC_DIMENSION]) # random seed used
    @test !isapprox(arr1, arr3)

    arr4 = as_array(wigner_state_sampled, [STOCHASTIC_DIMENSION, x]) # random seed used
    @test size(arr4) == (100, 64)

end


@testset "distinguishability of random variables" begin
    x = BasicDimension(RegularGrid(-5, 5, 32))
    rv1 = FluxionFields.random_normal([x], 0, 1/2)
    rv2 = FluxionFields.random_normal([x], 0, 1/2)

    arr1 = as_array(sample(rv1 - rv1, 10), [x, STOCHASTIC_DIMENSION])
    @test isapprox(arr1, zeros(arr1))

    arr2 = as_array(sample(rv1 - rv2, 10), [x, STOCHASTIC_DIMENSION])
    @test !isapprox(arr2, zeros(arr2))

    # Same RV only evaluate to the same random arrays only inside a single sample() call
    arr1 = as_array(sample(rv1, 10) - sample(rv1, 10), [x, STOCHASTIC_DIMENSION])
    @test !isapprox(arr1, zeros(arr1))

end


@testset "returned shape" begin
    x = BasicDimension(RegularGrid(-5, 5, 64))
    y = BasicDimension(RegularGrid(-5, 5, 32))

    @test size(as_array(x, [x, y])) == (64, 32)
    @test size(as_array(1, [x, y])) == (64, 32)

    dims = [x, y, STOCHASTIC_DIMENSION]
    @test size(as_array(sample(1, 10), dims)) == (64, 32, 10)
    @test size(as_array(sample(x, 10), dims)) == (64, 32, 10)
    @test size(as_array(sample(FluxionFields.random_normal([x], 0, 1), 10), dims)) == (64, 32, 10)
    @test size(as_array(sample(y + FluxionFields.random_normal([x], 0, 1), 10), dims)) == (64, 32, 10)
end
=#

#=
@testset "Wiener process" begin
    x = BasicDimension("x", RegularGrid(0, 1, 20))
    y = BasicDimension("y", RegularGrid(0, 1, 30))
    t = UnknownDimension("t")

    z = wiener_normal(x, t)

    func = as_array_function(x + y + z)
    func(x_arr, y_arr, dt, seed=123)

    # Need to handle:
    # - varying dt
    # - single/double step or some kind of noise splitting (check Num.Rec.)
    # - persistent state (for performance reasons)

end
=#

end
