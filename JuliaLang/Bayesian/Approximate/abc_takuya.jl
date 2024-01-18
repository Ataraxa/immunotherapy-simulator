using ApproxBayes
using DelimitedFiles: readdlm
using Distributions
using Pipe
using Random

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_restricted.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/treatments_lib.jl")
include("./distance.jl")
include("../../CommonLibrary/data_extractor.jl")
Random.seed!(1)

# Settings
input_treatment = CBD_IL_12_ver7
distro = Cauchy
data_set = 1

# Constants
dde_prob = create_problem(treatment=input_treatment)
prior_distro = @pipe distro |> Symbol(_) |> getfield(Main,_) 
fake_data_path = "Data/fakeData/trajectories-$data_set.csv"
days_array = input_treatment.active_days

# Simulator 
function tumour_growth_sim(params, constants, targetdata)
    p, u0 = repack_params(updateParams3(params...); do_split=true)
    pred = solve(dde_prob; p=p)
    pred_vol = @pipe pred[4,:] + pred[5,:] |> select_days(_, days_array, 0.1)

    return(distance(pred_vol, targetdata), 1)
end

# Target data
data_mat = readdlm(fake_data_path, ',')
data_mat = data_mat[:, days_array*trunc(Int, 1/0.1) .+ 1] # slice pred

# Inference 
setup = ABCSMC(
    tumour_growth_sim, # simulation function
    3,                 # number of parameters
    0.6                # target ϵ
    Prior([
        truncated(prior_distro(0.0, 1.0); lower=-100, upper=0),
        truncated(prior_distro(0.0, 1.0); lower=0   , upper=7),
        truncated(prior_distro(0.0, 1.0); lower=-100, upper=0)
    ])
)

smc = runabc(setup, targetdata; verbose=true, progress=true)