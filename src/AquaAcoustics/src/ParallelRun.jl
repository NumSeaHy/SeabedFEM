using Distributed
using DistributedArrays
using JLD2


function split_into_chunks(arr, num_chunks)
    chunk_size = ceil(Int, length(arr) / num_chunks)
    return [arr[i:min(i+chunk_size-1, end)] for i in 1:chunk_size:length(arr)]
end

function perform_distributed_computation(file_name)
    
    rmprocs(workers())
    addprocs(4)

    # Frequency range is fixed
    N = 100 # Number of steps
    freq_range_start = 1.0e3
    freq_range_stop = 20.0e3
    freq_range_step = (freq_range_stop-freq_range_start) / (N-1)

    # Define and distribute the frequency range
    freq_range = freq_range_start:freq_range_step:freq_range_stop


    # Include the function definition on all workers
    @everywhere include("./RunParallelFunc.jl")
    @everywhere filename = $file_name  # Ensure filename is available on all workers

    # Fetch local parts of the distributed array
    freq_chunks = split_into_chunks(freq_range, nworkers())

    # Prepare arguments for pmap, passing filename with each frequency chunk
    args = [(chunk, filename) for chunk in freq_chunks]

    # Execute Run function in parallel
    results = vcat(pmap(x -> Run(x[1], x[2]), args)...)

    # Save the results
    @save "results/"*file_name*".jld2" freq_range results

    return results
end


filename = "scenario1"
perform_distributed_computation(filename)
rmprocs(workers())
# Function to perform the distributed computation
