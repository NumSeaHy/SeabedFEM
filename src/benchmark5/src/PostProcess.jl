"""
This file opens the jld2 files with the solutions of the two simulations and compares them 
computing the L2 norm of the error. It also writes the results in a .vtu file that can be
open in an unique Paraview session.
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using Gridap.Io
using Gridap.CellData
# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
include("./AnalyticalIncident.jl")

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Define the tags of the mesh boundary
dimension = 2
labels = get_face_labeling(model)
tags = get_face_tag(labels, dimension)

# Load the Single.FE.func
Uh1 = from_jld2_file("./results/FemSol.jld2") # Conventional
Uh2 = from_jld2_file("./results/FemSolTrans.jld2") # Translated

# Get the tiangulation and the physical coordinates of the mesh
Ω = get_triangulation(model)
xp = get_physical_coordinate(Ω)
Ω_Ph = Triangulation(model, tags="physical_domain")
degree = 2
dΩ_Ph = Measure(Ω_Ph, degree) 

# Take the DOF values of the Single.FE.func == solution
val1 = get_cell_dof_values(Uh1)
val2 = get_cell_dof_values(Uh2)

# Fucntions to make the post-processing of the solution
γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
γ_2(x) = x[2] > 0 ? 1. : 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 
H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0.0+0.0im, 0.0+0.0im, γ_2(x))
J(x) = det(H(x))
Jinv(x) = 1/det(H(x))

# PML boundary condition incident field
# A_F, B_F, A_P, B_P, A_PML, B_PML = solve_coefficients_pml(ω, ρ_F(ω)*c_F(ω), ρ_P(ω)*c_P(ω), P_0, t_P, t_F, k_F(ω), k_P(ω), σ_0, d_PML)
# u_incident(x) = exact_quadratic_pml(x, t_P, k_F(ω), k_P(ω), d_PML, σ_0, A_F, B_F, A_P, B_P, A_PML, B_PML)
# Sommerfeld boundary condtion incident field
u_incident(x) = exact_solution_som(x, ρ_F(ω), c_F(ω), k_F(ω), ρ_P(ω), c_P(ω), k_P(ω), P_0, t_F, t_P, d_PML, σ_0)

# Define the finite element space previously used in perform the two simulations: Raviart-Thomas of order 1 with the corresponding boundary conditions
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom", "rock"], vector_type=Vector{ComplexF64})

# Generate the CellFields for the two solutions using the Finite Element Space and the DOF values
s1 = CellField(V, val1) # Conventional
s2 = CellField(V, val2) # Translated

# Post-Processing the solution
s1 = (Jinv ∘ xp) * ((H ∘ xp) ⋅ s1)
u_i = u_incident ∘ xp
s2 = (Jinv ∘ xp) * ((H ∘ xp) ⋅ s2) + u_i

# Not neccesarly to compute the L2 norm of the error, but I already had the code from another example
s1x = (u->u[1]) ∘ s1
s1y = (u->u[2]) ∘ s1
s2x = (u->u[1]) ∘ s2
s2y = (u->u[2]) ∘ s2

e = (s2-s1)./s1

# Write the solution in a .vtu file
# writevtk(Ω_Ph,"./results/result_compare.vtu", cellfields=["Re(uh1)"=>real(s1), "Im(uh1)"=>imag(s1),
                                                        #   "Re(uh2)"=>real(s2), "Im(uh2)"=>imag(s2),
                                                        #   "Re(err)"=>real(e),  "Im(err)"=>imag(e)])     
# Compute the L2 norm of the relative error
100 * sqrt(sum(∫(abs2((s2x-s1x)))*dΩ_Ph)/sum(∫(abs2((s2x)))*dΩ_Ph)+
           sum(∫(abs2((s2y-s1y)))*dΩ_Ph)/sum(∫(abs2((s2y)))*dΩ_Ph))



