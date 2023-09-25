using Evolutionary
using JLD2

res = Evolutionary.optimize(x->-sum(x),
                        BitVector(zeros(30)),
                        GA(selection=uniformranking(5),
                        mutation=flip, crossover=SPX),
                        Evolutionary.Options(iterations=10))

save_object("dummy_goug.jld2", res)