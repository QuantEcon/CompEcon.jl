module CompEconTests

using CompEcon

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

tests = [
    "types.jl"
    ]

if length(ARGS) > 0
    tests = ARGS
end

end_jl(s) = endswith(s, ".jl") ? s : s * ".jl"

for t in tests
    print_with_color(:green, "* $t\n")
    include(end_jl(t))
end


end  # module
