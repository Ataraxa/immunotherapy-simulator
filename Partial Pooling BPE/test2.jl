using Pkg
Pkg.status()
Pkg.add("StatsBase")

p1 = [1, 2, 3]
p2 = [4, 5, 6]

p = [p1, p2]
print(p[1][2])