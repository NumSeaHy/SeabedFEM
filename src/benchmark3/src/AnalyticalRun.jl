using Gridap

include("./Configuration.jl")
include("./AnalyticalSommerfeld.jl")
include("NumericalSystemSolve.jl")

model = DiscreteModelFromFile("data/pml_mesh.json")

A_F, B_F, A_P, B_P, A_PML, B_PML = solve_coefficients_pml(ω, Z_F, Z_P, P_0, t_P, t_F, k_F, k_P, σ_3, d_PML)

Ω = Triangulation(model)

u1 = CellField(x-> ϕ_pml(x, t_P, k_F, k_P, σ_3, A_F, B_F, A_P, B_P, A_PML, B_PML), Ω)


writevtk(Ω,"pml_analytical.vtu", cellfields=["Re_y"=>real(u1),
                                             "Im_y"=>imag(u1)])

# writevtk(Ω,"pml_analytical.vtu", cellfields=["RealSom_y"=>real(u1),
#                                              "ImagSom_y"=>imag(u1),
#                                              "RealNorm_y"=>real(u2),
#                                              "ImagNorm_y"=>imag(u2)])
