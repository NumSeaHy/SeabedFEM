using CairoMakie
using JLD2
function Plot(path, savename)
    @load path errors logspace 
    f = Figure(size = (500,200), fontsize = 13)
    ax = Axis(f[1, 1],
    xlabel = L"\sigma",
    ylabel = L"100\frac{‖{u_h-u}‖_{L²(Ω)}}{‖{u}‖_{L²(Ω)}}",
    xscale = log10,
    )
    lines!(ax, logspace, errors,  color = :black)
    # scatter!(ax, logspace, errors, markersize=10, marker = :circle)
    save("./results/"*savename*".pdf", f)
end