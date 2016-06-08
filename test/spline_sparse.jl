# Tests constructors and subtype relationships
module TestSplineSparse

using CompEcon
using Base.Test
using FactCheck

facts("Test SplineSparse") do

    context("full") do
        s = SplineSparse{1,2,Int,Int}(3, [1, 2], [1, 2, 3, 4])
        want = [1 2 0; 0 3 4]
        @fact full(s) --> want

        s = SplineSparse{2,2,Int,NTuple{2,Int}}(6, [(1, 5), (2, 4)], 1:8)
        want = [1 2 0 0 3 4
                0 5 6 7 8 0]

        @fact full(s) --> want
    end

    context("row_kron") do
        s1 = SplineSparse{1,2,Int,Int}(3, rand(1:2, 6), rand(1:10, 12))
        s2 = SplineSparse{1,2,Int,Int}(4, rand(1:3, 6), rand(1:10, 12))
        want = row_kron(full(s1), full(s2))
        s12 = row_kron(s1, s2)

        @fact full(s12) --> want
        @fact full(row_kron(s12, s1)) --> row_kron(full(s12), full(s1))
        @fact full(row_kron(s12, s12)) --> row_kron(full(s12), full(s12))

        Base.zero(::Type{ASCIIString}) = ""

        s1 = SplineSparse{1,2,ASCIIString,Int}(2, [1, 1], ["a", "b", "c", "d"])
        s2 = SplineSparse{1,3,ASCIIString,Int}(3, [1, 1], map(string, 1:6))

        want = ["a1" "a2" "a3" "b1" "b2" "b3"
                "c4" "c5" "c6" "d4" "d5" "d6"]
        have = row_kron(s1, s2)

        @fact full(have) --> want

    end

    context("getindex") do
        s = SplineSparse{1,2,Float64,Int}(10, rand(1:9, 6), rand(12))
        full_s = full(s)

        for r in 1:size(s, 1)
            for c in 1:size(s, 1)
                @fact s[r, c] --> full_s[r, c]
            end
        end
    end

    context("*") do
        s1 = SplineSparse{1,2,Int,Int}(3, [1, 2], 1:4)
        s2 = SplineSparse{1,3,Int,Int}(4, [2, 1], 1:6)
        s12 = SplineSparse{2,3,Int,NTuple{2,Int}}(12, [(2, 6), (5, 9)],
                            [1, 2, 3, 2, 4, 6, 12, 15, 18, 16, 20, 24])

        x1 = rand(3)
        x2 = rand(4)
        x12 = rand(12)

        @fact s1*x1 --> full(s1)*x1
        @fact s2*x2 --> full(s2)*x2
        @fact s12*x12 --> full(s12)*x12
    end

end

end  # module
