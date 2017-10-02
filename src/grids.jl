abstract type Grid end


function _is_monotonic_increasing{T}(x::Array{T, 1})
    if length(x) < 2
        return true
    end

    all(elem > 0 for elem in x[2:end] .- x[1:end-1])
end


struct RegularGrid <: Grid
    left :: Float64
    right :: Float64
    points :: Int
    vertices :: Array{Float64, 1}

    function RegularGrid(left=0, right=1, points=50)
        return new(left, right, points, linspace(left, right, points))
    end
end


struct ArbitraryGrid <: Grid
    vertices :: Array{Float64, 1}
    monotonic :: Bool

    ArbitraryGrid(data) = new(data, _is_monotonic_increasing(data))
end


grid_step(g::RegularGrid) = (g.right - g.left) / g.points


vertices(g::RegularGrid) = g.vertices
vertices(g::ArbitraryGrid) = g.vertices


is_monotonic(g::RegularGrid) = true
is_monotonic(g::ArbitraryGrid) = g.monotonic
