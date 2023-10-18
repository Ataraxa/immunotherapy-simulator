struct Test 
    a
    b
end
test = Test(1, 2)

mutable struct basedStruct 
    a
    b
    c
    d
end
based = basedStruct(10, 20, 30, 40)

for name in fieldnames(typeof(test))
    field = Symbol(name)
    code = quote
        (obj) -> obj.$field 
    end
    val = eval(code)(test)
    setproperty!(based, field, val)
end