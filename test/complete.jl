module TestComplete

using CompEcon
using FactCheck

facts("Testing logic using strings") do
    Base.one{T<:AbstractString}(::Type{T}) = "one "
    z = ["x " "y " "z "]

    want_1 = ["one " "x " "y " "z "]
    @fact complete_polynomial(z, 1) --> want_1

    want_2 = ["one " "x " "y " "z " "one x x " "one x y " "one x z " "one y y " "one y z " "one z z "]
    @fact complete_polynomial(z, 2) --> want_2
end

facts("testing numerical routintes") do
    z = rand(100, 10)
    the_ones = ones(size(z, 1))
    for d in 1:5
        n_comp = n_complete(size(z, 2), d)
        buff_d = Array(Float64, size(z, 1), n_comp)
        out_d = complete_polynomial(z, d)
        complete_polynomial!(z, d, buff_d)
        @fact size(out_d, 1) --> size(z, 1)
        @fact size(out_d, 2) --> n_comp
        @fact out_d[:, 1] --> the_ones
        @fact out_d[:, 2:size(z, 2)+1] --> z
        @fact out_d --> buff_d
    end

    z = rand(10, 3)
    out_2 = complete_polynomial(z, 2)
    # need to test columns size(z, 2) + 2:end
    @fact out_2[:, 5] --> z[:, 1].^2
    @fact out_2[:, 6] --> z[:, 1].*z[:, 2]
    @fact out_2[:, 7] --> z[:, 1].*z[:, 3]
    @fact out_2[:, 8] --> z[:, 2].*z[:, 2]
    @fact out_2[:, 9] --> z[:, 2].*z[:, 3]
    @fact out_2[:, 10] --> z[:, 3].*z[:, 3]
end

end
