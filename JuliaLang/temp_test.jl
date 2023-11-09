a = [1 0 0; 1 0 0; 0 0 0]

a = a[:, vec(mapslices(col -> any(col .!= 0), a, dims = 1))] 
a = a[vec(mapslices(col -> any(col .!= 0), a, dims = 2)), :] 