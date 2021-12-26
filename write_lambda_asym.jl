using JLD2, DelimitedFiles
include("/scratch/projects/hhp00048/codes/lDGAPostprocessing/src/IO.jl")

fname = "lambda_asym.jld2"

jldopen(fname, "w") do f
    a,b,c = read_λ_asympt("trilex_omega")
    f["γup"] = a; f["γdo"] = b; f["γpp"] = c;
    a,b,c = read_λ_asympt("tripamp_omega")
    f["λup_amp"] = a; f["λdo_amp"] = b; f["λpp_amp"] = c;
    a,b,c = read_λ_asympt("trip_omega")
    f["λup"] = a; f["λdo"] = b; f["λpp"] = c;
end;
