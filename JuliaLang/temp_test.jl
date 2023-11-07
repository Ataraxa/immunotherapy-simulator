using Match

a = "double"
b = 10
@match a begin
    "double" => begin
        b = 2
        println("haha")
    end
    "triple" => (b=3)
end

println(b)