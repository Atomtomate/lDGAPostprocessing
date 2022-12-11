using JLD2, DelimitedFiles
include("/scratch/projects/hhp00048/codes/lDGAPostprocessing/src/IO.jl")

fname = "lambda_asym.jld2"

jldopen(fname, "w") do f
    a,b,c = read_λ_asympt("trip_omega")
    f["λch"] = a .+ b; f["λsp"] = a .- b; f["λpp"] = c;
    a,b,c = read_λ_asympt("trilex_omega")
    f["γch"] = a; f["γsp"] = b; f["γpp"] = c;
    a,b,c = read_λ_asympt("tripamp_omega")
    f["λch_amp"] = a; f["λsp_amp"] = b; f["λpp_amp"] = c;
end;
