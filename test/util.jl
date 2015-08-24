module TestUtil

using CompEcon
using FactCheck

facts("test CompEcon.lookup") do
    table1 = [1.0, 4.0]
    table2 = [1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0]

    x = [0.5, 1.0, 1.5, 4.0, 5.5]
    x2 = [0.5, 2.0]
    @fact CompEcon.lookup(table1, x, 0) --> [0, 1, 1, 2, 2]
    @fact CompEcon.lookup(table1, x, 1) --> [1, 1, 1, 2, 2]
    @fact CompEcon.lookup(table1, x, 2) --> [0, 1, 1, 1, 1]
    @fact CompEcon.lookup(table1, x, 3) --> [1, 1, 1, 1, 1]

    @fact CompEcon.lookup(table2, x, 0) --> [0, 3, 3, 12, 12]
    @fact CompEcon.lookup(table2, x, 1) --> [3, 3, 3, 12, 12]
    @fact CompEcon.lookup(table2, x, 2) --> [0, 3, 3, 8, 8]
    @fact CompEcon.lookup(table2, x, 3) --> [3, 3, 3, 8, 8]


    @fact CompEcon.lookup([1.0], x2, 0) --> [0, 1]
    @fact CompEcon.lookup([1.0], x2, 1) --> [1, 1]
    @fact CompEcon.lookup([1.0], x2, 2) --> [0, 1]
    @fact CompEcon.lookup([1.0], x2, 3) --> [1, 1]

    # test scalar version of lookup
    x2 = collect(linspace(-2.0, 4.0, 10))
    @fact [CompEcon.lookup(x2, -3.0, i) for i=0:3] --> [0, 1, 0, 1]
    @fact [CompEcon.lookup(x2, 5.0, i) for i=0:3] --> [10, 10, 9, 9]
    @fact [CompEcon.lookup(x2, i, 0) for i=x2] --> collect(0:length(x2)-1)
    @fact [CompEcon.lookup(x2, i, 1) for i=x2] --> [1; 1:length(x2)-1]
end

facts("test CompEcon._check_order") do
    # check (::Int, ::Int) method
    @fact CompEcon._check_order(10, 0) --> fill(0, 1, 10)
    @fact CompEcon._check_order(1, 0) --> fill(0, 1, 1)

    # check (::Int, ::Vector) method
    ov = [0, 0]
    @fact CompEcon._check_order(2, ov) --> reshape(ov, 1, 2)
    @fact_throws DimensionMismatch CompEcon._check_order(3, ov)

    # check (::Int, ::Matrix) method
    om = [0 0]
    @fact CompEcon._check_order(2, om) --> om
    @fact CompEcon._check_order(1, om) --> om'
    @fact_throws DimensionMismatch CompEcon._check_order(3, ov)
end


facts("test CompEcon.gridmake") do
    # TODO: these tests are very incomplete. They just test that gridmake
    #       preserves input type

    for T in [Float16, Float32, Float64, Int128, Int16, Int32, Int64, Int8]
        @fact eltype(CompEcon.gridmake(rand(T, 2), rand(T, 2))) --> T
    end
end


facts("test CompEcon.ckronx") do
    # will test by constructing an interpoland, then evaluating at the nodes
    # and verifying that we get back our original function
    basis = CompEcon.Basis(CompEcon.Basis(CompEcon.Spline(), [-1.0, 1.0], 13, 3),
                           CompEcon.Basis(CompEcon.Spline(), [-5.0, 3.0], 18, 3))
    X, x12 = CompEcon.nodes(basis);

    # make up a funciton and evaluate at the nodes
    f(x1, x2) = cos(x1) ./ exp(x2)
    f(X::Matrix) = f(X[:, 1], X[:, 2])
    y = f(X)

    # fit the interpoland in Tensor form (tensor b/c using x12)
    c, bs = CompEcon.funfitxy(basis, x12, y);

    # verify that we are actually interpolating -- all heavy lifting in funeval
    # is done by ckronx so this is effectively testing that we wrote that
    # function properly
    @fact funeval(c, bs, [0 0]) --> roughly(y; atol=1e-14)
end

end  # module
