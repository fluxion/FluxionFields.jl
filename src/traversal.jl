type TraversalControls
    skip_children :: Bool
    call_after :: Bool
    called_after :: Bool
    args :: Any
    TraversalControls(args) = new(false, false, false, args)
end


function skip_children!(controls)
    controls.skip_children = true
end


function call_after!(controls)
    controls.call_after = true
end


function replace_args!(controls, args...)
    controls.args = args
end


function inspect_and_modify(func::Function, state, val, args...)

    controls = TraversalControls(args)

    state, new_val = func(controls, state, val, args...)

    children_args = controls.args

    if !controls.skip_children && new_val === val
        subnodes = to_subnodes(val)
        if !(subnodes === nothing)
            new_subnodes = Array{Any, 1}(length(subnodes))
            changed = false
            for (i, subnode) in enumerate(subnodes)
                state, new_subnode = inspect_and_modify(func, state, subnode, children_args...)
                changed |= !(new_subnode === subnode)
                new_subnodes[i] = new_subnode
            end

            new_val = changed ? from_subnodes(val, new_subnodes) : val
        end
    end

    if controls.call_after
        controls.called_after = true
        state, new_val = func(controls, state, new_val, args...)
    end

    state, new_val
end


function inspect(func::Function, state, val, args...)
    state, _ = inspect_and_modify(state, val, args...) do controls, state, val, args...
        func(controls, state, val, args...), val
    end
    state
end


function modify(func::Function, val, args...)
    _, new_val = inspect_and_modify(nothing, val, args...) do controls, state, val, args...
        state, func(controls, val, args...)
    end
    new_val
end
