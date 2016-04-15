test_files = ["types.jl", "basis.jl", "optimization.jl", "util.jl",
              "spline.jl",	"interp.jl", "cheb.jl", "lin.jl",
              "basis_structure.jl", "zeros.jl", "complete.jl"]

for f in test_files
    include(f)
end
