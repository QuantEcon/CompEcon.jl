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

end


facts("test CompEcon.gridmake") do
    # TODO: these tests are very incomplete. They just test that gridmake
    #       preserves input type

    for T in [Float16, Float32, Float64, Int128, Int16, Int32, Int64, Int8]
        @fact eltype(CompEcon.gridmake(rand(T, 2), rand(T, 2))) --> T
    end
end

end  # module
