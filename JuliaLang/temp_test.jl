using Pipe

arr = @pipe zeros(3) |> exp.(_[2])