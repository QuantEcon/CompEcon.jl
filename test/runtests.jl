test_files = ["types.jl", "basis.jl", "optimization.jl", "test_util.jl"]

for f in test_files
    include(f)
end
