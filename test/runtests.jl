test_files = ["types.jl", "basis.jl", "optimization.jl", "util.jl",
              "quad.jl",
              "basis_structure.jl"]

for f in test_files
    include(f)
end
