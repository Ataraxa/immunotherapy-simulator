# using Match

# function prob_fact(model)
#     a::Function = (x->print(x))      
#     @match model begin
#         "double" => begin
#             a = 2
#         end

#         "triple" => begin 
#             a = function(input)
#                 return input * 3
#             end
#         end
#     end

#     return a
# end

# a = prob_fact("triple")
# a(3)
using ForwardDiff: Dual

p = [
    Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(17.493801636164296, 0.0, 2.339124341340631, 0.0),
    Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(0.2785968912154963, 0.0, 0.0, 1.141240725012459)
]
new_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p)
        println(typeof(param))
        if typeof(param) <: Dual
            new_p[i] = param.value
        end
    end
new_p

