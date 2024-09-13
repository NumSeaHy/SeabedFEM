using CairoMakie
using JLD2

function PlotFrequencyResp(path)
    @load path f L2
    fig = Figure(size = (500, 200), fontsize = 15)
    ax = Axis(fig[1, 1],
    xlabel = L"f \, [\mathrm{kHz}]",
    ylabel = L"\sqrt{∑_{i=1}^{2}∫|u_i|^2dΩ}")
    lines!(ax, f/1e3, L2, color=:black)
    # scatter!(ax, f, L2, markersize=10, marker = :circle, color=:black)
    save("./results_freq/freq_plot.pdf", fig)
    fig    
end


function vectors_to_txt(file_path::String, vectors::Vector{Vector{T}}) where T
    # Abrir el archivo para escritura
    open(file_path, "w") do file
        # Compute the maximum length of the vectors
        max_length = maximum(length.(vectors))
        
        # Iterate by each index to obtain the elements of the vectors
        for i in 1:max_length
            # Para cada vector, obtener el elemento si el índice es válido, o usar un espacio en blanco si no lo es
            line = [i <= length(vec) ? string(vec[i]) : "" for vec in vectors]
            
            # Join the elemetns and write the file
            write(file, join(line, " ") * "\n")
        end
    end
end

