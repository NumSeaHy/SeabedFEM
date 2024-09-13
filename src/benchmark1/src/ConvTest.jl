using CairoMakie
include("Run.jl")

h = [0.06, 0.06/2, 0.06/4, 0.06/8]
el2 = Float64[]
meshes = ["data/benchmark_2_coarse", "data/benchmark2_medium", "./data/benchmark2_fine", "./data/benchmark2_extremely_fine"]

for i in meshes
    push!(el2, run(i))
end

x = log10.(h)
y = log10.(el2)
linreg = hcat(fill!(similar(x), 1), x) \ y
m = linreg[2]
println("Convergence order is: $m")

f = Figure(size = (600,300), fontsize = 11)
ax = Axis(f[1, 1],
    xlabel = L"h",
    ylabel = L"100\frac{‖{u_h-u}‖_{L²(Ω)}}{‖{u}‖_{L²(Ω)}}",
    xscale = log10,
    yscale = log10
)
lines!(ax, h, el2)
scatter!(ax, h, el2, markersize=10, marker = :circle)
save("./results/order_bm2.pdf", f)