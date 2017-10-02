import Base: show


function show(io::IO, x::UnknownDimension)
    print(io, "UD[$(x.id)]")
end

function show(io::IO, x::BasicDimension)
    print(io, "BD[$(x.id)](")
    show(io, x.grid)
    print(io, ")")
end

function show(io::IO, x::RegularGrid)
    print(io, "UG(")
    show(io, x.left)
    print(io, ", ")
    show(io, x.right)
    print(io, ", ")
    show(io, x.points)
    print(io, ")")
end

function show(io::IO, x::TransformedDimension)
    print(io, "TD[$(x.id)](")
    show(io, x.base_dimension)
    print(io, ")")
end

function show(io::IO, x::UnknownField)
    print(io, "UF[$(x.id)](")
    for (i, dim) in enumerate(x.dimensions)
        show(io, dim)
        if i != length(x.dimensions)
            print(io, ", ")
        end
    end
    print(io, ")")
end

function show(io::IO, x::KnownField)
    print(io, "KF(")
    if length(x.dimensions) > 0
        for (i, dim) in enumerate(x.dimensions)
            show(io, dim)
            if i != length(x.dimensions)
                print(io, ", ")
            end
        end
    else
        show(io, reshape(from_backend(x.data), 1)[1])
    end
    print(io, ")")
end

function show(io::IO, x::Fourier)
    if x.inverse
        print(io, "I")
    end
    print(io, "Fourier[")
    for (i, (from, to)) in enumerate(x.dimensions)
        show(io, from)
        print(io, "->")
        show(io, to)
        if i != length(x.dimensions)
            print(io, ", ")
        end
    end
    print(io, "](")
    show(io, x.fexpr)
    print(io, ")")
end


function show(io::IO, fexpr::FieldExprNode)
    print(io, "{")
    show(io, fexpr.head)
    print(io, ", ")

    args = fexpr.args
    for (i, arg) in enumerate(args)
        show(io, arg)
        if i != length(args)
            print(io, ", ")
        end
    end
    print(io, "}")
end
