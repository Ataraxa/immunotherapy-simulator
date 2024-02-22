using Combinatorics

x = collect(1:9)
comx = combinations(x, 2)
for c in comx
    println(c)
end