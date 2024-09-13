using Gridap
using Gridap.Fields
using Gridap.Geometry

include("./Configuration.jl")
include("./NumericalSystemSolve.jl")

model = DiscreteModelFromFile("data/pml_mesh.json")

labels = get_face_labeling(model)
dimension = 2
tags = get_face_tag(labels, dimension)

fluid_PML_tag = get_tag_from_name(labels,"fluid_PML")
porous_PML_tag = get_tag_from_name(labels,"porous_PML")
bottom_PML_tag = get_tag_from_name(labels,"bottom_PML")
bottom_PML_lat_tag = get_tag_from_name(labels,"bottom_PML_lat") # Corners PML 
fluid_tag = get_tag_from_name(labels,"fluid_domain")
porous_tag = get_tag_from_name(labels,"porous_domain")
transducer_tag = get_tag_from_name(labels,"transducer")
transducer_pml_tag = get_tag_from_name(labels,"transducer_pml")  

order = 1
V = TestFESpace(model, ReferenceFE(raviart_thomas, Float64, order),
    conformity=:HDiv, dirichlet_tags=["sides"], vector_type=Vector{ComplexF64})

uD = VectorValue(0.0, 0.0)

U = TrialFESpace(V, [uD])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

k_F = 2 * pi * f/c_F
k_P = 2 * pi * f/c_P 

neumanntags = ["transducer", "transducer_pml"]
Γ = BoundaryTriangulation(model, tags="transducer")
dΓ = Measure(Γ, degree)
Γ_pml = BoundaryTriangulation(model, tags="transducer_pml")
dΓ_pml = Measure(Γ_pml, degree)
nb = get_normal_vector(Γ) # Normal vector to the transducer boundary
nb2 = get_normal_vector(Γ_pml) # Normal vector to the transducer boundary
τ = CellField(tags, Ω) # Tags field
p = P_0 # Transducer pressure
σ_1F = σ_3
σ_1P = σ_3
σ_2P = σ_3
γ_1_F = 1 + 1im/k_F * (σ_1F)
γ_1_P = 1 + 1im/k_P * (σ_1P)
γ_2_P = 1 + 1im/k_P * (σ_2P) # +

function ρ(tag)
      if tag == fluid_tag || tag == fluid_PML_tag
            return ρ_F
      else
            return ρ_P
      end
end

function K(tag)
      if tag == fluid_tag || tag == fluid_PML_tag
            return ρ_F*c_F^2
      else
            return ρ_P*c_P^2
      end
end

function pml_stiffness(tag, ∇u, ∇v)
      # return (TensorValue{2,2,ComplexF64}(\gamma_2(tag), 0.0+0.0im, 0.0+0.0im, \gamma_1) ⋅ (TensorValue{2,2,ComplexF64}(1/\gamma_1, 0.0+0.0im, 0.0+0.0im, 1/\gamma_2) ⊙ ∇u))     
      if tag == fluid_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_F, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⊙ ∇u) ⋅ TensorValue{2,2, ComplexF64}(1.0+0.0im, 0.0+0.0im, 0.0+0.0im, γ_1_F) ⊙ ∇v)
      elseif tag == porous_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_P, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⊙ ∇u) ⋅ (TensorValue{2,2, ComplexF64}(1.0+0.0im, 0.0+0.0im, 0.0+0.0im, γ_1_P) ⊙ ∇v))
      elseif tag == bottom_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1.0+0.0im, 0.0+0.0im, 0.0+0.0im, 1/γ_2_P) ⊙ ∇u) ⋅ (TensorValue{2,2, ComplexF64}(γ_2_P, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⊙ ∇v))
      elseif tag == bottom_PML_lat_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_P, 0.0+0.0im, 0.0+0.0im, 1/γ_2_P) ⊙ ∇u) ⋅ (TensorValue{2,2,ComplexF64}(γ_2_P, 0.0+0.0im, 0.0+0.0im, γ_1_P) ⊙ ∇v))
      else
            return ((TensorValue{2,2,ComplexF64}(1.0+0.0im, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⊙ ∇u) ⋅ (TensorValue{2,2,ComplexF64}(1.0+0.0im, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⊙ ∇v))
      end
end

function pml_mass(tag)
      if tag == fluid_PML_tag
            return γ_1_F
      elseif tag == porous_PML_tag
            return γ_1_P
      elseif tag == bottom_PML_tag
            return γ_2_P
      elseif tag == bottom_PML_lat_tag
            return γ_1_P * γ_2_P
      else
            return 1.0 + 0.0im
      end
end

# Bilinear term
a(u, v) = ∫( (ρ ∘ (τ)) * (pml_mass ∘ τ) * (u⋅v) * ω^2 )dΩ - ∫( (K ∘ (τ)) * (pml_stiffness ∘ (τ, gradient(u), gradient(v))) )dΩ
# Linear temr
#b(v) = ∫( (pml_source ∘ (τ, v, nb)) * p)dΓ
 
# b(v) = ∫( (v⋅nb) * p)dΓ + ∫( ((TensorValue{2,2,ComplexF64}(1/γ_1_F, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⋅ v)⋅nb) * p)dΓ_pml
 
b(v) = ∫( (v⋅nb) * p)dΓ + ∫( ((TensorValue{2,2,ComplexF64}(1/γ_1_F, 0.0+0.0im, 0.0+0.0im, 1.0+0.0im) ⋅ v)⋅nb2) * p)dΓ_pml
# b(v) = ∫( (v⋅nb) * p)dΓ
# Assembly the system
op = AffineFEOperator(a, b, U, V)
# Solve the system
uh = solve(op)
# Post-Procesing
uh_x = (u->u[1])∘uh
uh_y = (u->u[2])∘uh
function AnalyticalBox(tag) # Function to only quantifies the error in the bottom PML
      if tag == bottom_PML_tag || tag == bottom_PML_lat_tag
      return 1.0
      else
      return 0.0
      end    
end
# Generate a field with the tags
λ = AnalyticalBox ∘ τ
# Compute the analytical solution
A_F, B_F, A_P, B_P, A_PML, B_PML = solve_coefficients_pml(ω, Z_F, Z_P, P_0, t_P, t_F, k_F, k_P, σ_3, d_PML)
u = CellField(x-> ϕ_pml(x, t_P, k_F, k_P, σ_3, A_F, B_F, A_P, B_P, A_PML, B_PML), Ω)
# Compute the error and apply the mask to the error
e = u - uh_y
e = λ ⋅ e

writevtk(Ω,"./results/bottom_pml_e.vtu", cellfields=[ "Re_y"=>real(uh_y),
                                                "Im_y"=>imag(uh_y),
                                                "Re_x"=>real(uh_x),
                                                "Im_x"=>imag(uh_x),
                                                "Im_analytical"=>imag(u),
                                                "Re_analytical"=>real(u),
                                                "Im_error"=>imag(e),
                                                "Re_error"=>real(e)])