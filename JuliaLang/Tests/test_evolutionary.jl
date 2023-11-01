using Evolutionary

res = Evolutionary.optimize(x->sum(x), BitVector(zeros(3)), GA(),
Evolutionary.Options(
    iterations = 20
))