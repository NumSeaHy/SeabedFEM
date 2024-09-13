"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl. The solution is saved in a .vtu file and
the Singled.FE.Function is saved in a .jld2 file to the compare  the results with the 
translated approach.
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Define the tags of the mesh boundary
dimension = 2
labels = get_face_labeling(model)
tags = get_face_tag(labels, dimension)

# Define the finite element space: Raviart-Thomas of order 1
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom", "rock"], vector_type=Vector{ComplexF64})

# Define the trial function with null Dirichlet boundary conditions
uD_null = VectorValue(0.0, 0.0)
U = TrialFESpace(V, [uD_null, uD_null, uD_null])

# Triangulation and spatial coordinates of the computational domain
Ω = Triangulation(model)
xp = get_physical_coordinate(Ω)

# Define the measure for the fluid and porous domains
degree = 2
dΩ = Measure(Ω, degree)

# Define the measure for the transducer boundary
Γ = BoundaryTriangulation(model, tags="transducer")
dΓ = Measure(Γ, degree)
n = get_normal_vector(Γ) # Normal vector to the transducer boundary

# Define the physical coefficients in the fluid and porous media: bulk modulus and mass density
K(x) = x[2]-t_P > 0 ? K_F(ω) : K_P(ω)
ρ(x) = x[2]-t_P > 0 ? ρ_F(ω) : ρ_P(ω) 

# Definition of the PML tensors H, H^-1, and the Jacobian for the fluid and porous PML (quadratic profile)
γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
γ_2(x) = x[2] > 0 ? 1. : 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 
H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
Hinv(x) = inv(H(x))
J(x) = det(H(x))
Jinv(x) = 1/det(H(x))

# Bilinear term
a(u, v) = ∫( (K ∘ xp) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
      ∫( (ρ ∘ xp) * (Jinv ∘ xp) * (((H ∘ xp) ⋅ u) ⋅ ((H ∘ xp) ⋅ v)) * ω^2 )*dΩ

# Source term
b(v) = ∫(-P_0 * (v ⋅ n))*dΓ

# Assembly the FE matrices and solve the linear system
op = AffineFEOperator(a, b, U, V)
Uh = solve(op)

# Post-Processing the solution in the PML layers (uh=Uh in the physical domain)
uh = (Jinv ∘ xp) * ((H ∘ xp) ⋅ Uh)

# Save the results in a vtk file and save the Singled.FE.Function in a .jld2 file
# writevtk(Ω,"./results/result.vtu", cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh)]) 
to_jld2_file(Uh , "./results/FemSol.jld2") 


