using Match

a = "double"
@match a begin
    "double" => global b = sqrt(6)
    "triple" => (b=3)
end

println(b)