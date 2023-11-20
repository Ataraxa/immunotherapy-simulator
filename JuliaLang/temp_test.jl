using Match 

cdt = "logarithm"

@match cdt begin 
    "identity" => global transform = (x -> x)
    "logarithm" => global transform = (x -> log.(x))
end

transform([1, 2, 3])