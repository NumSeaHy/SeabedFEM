using Peaks
using JLD2
include("FrequencyResponse.jl")

@load "results_freq/freq.jld2" f L2

n = 3
pks, vals = findmaxima(L2)
aux = zeros(n)
for i in 1:n
    aux[i] = f[pks[i]]
end

Run(aux, true, false)