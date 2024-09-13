using Gridap
using Gridap.Geometry
using Gridap.Fields
using Gridap.Io
using GridapGmsh

include("./Properties.jl")
include("./Configuration.jl")
include("./Analytical.jl")

function run(mesh)
      
      function ρ(tag)
            if tag == water_tag
                  return ρ_F
            else
                  return ρ_P
            end
      end

      function K(tag)
            if tag == water_tag
                  return ρ_F*c_F^2
            else
                  return ρ_P*c_P^2
            end
      end
      
      
      model = DiscreteModelFromFile(mesh*".json")

      labels = get_face_labeling(model)
      dimension = 2
      tags = get_face_tag(labels, dimension)

      water_tag = get_tag_from_name(labels,"water_domain")
      order = 1
      V = TestFESpace(model, ReferenceFE(raviart_thomas, Float64, order),
            conformity=:HDiv, dirichlet_tags=["left_water","right_water","left_porous",
                                              "right_porous", "bottom"],vector_type=Vector{ComplexF64})

      uD = VectorValue(0.0, 0.0)
      U = TrialFESpace(V, [uD, uD, uD, uD, uD])
      degree = 2
      Ω = Triangulation(model)
      dΩ = Measure(Ω, degree)

      neumanntags = ["transducer"]
      Γ = BoundaryTriangulation(model, tags=neumanntags)
      dΓ = Measure(Γ, degree)

      nb = get_normal_vector(Γ)
      
      τ = CellField(tags, Ω)

      p = P_0
      a(u, v) = ∫((ρ∘(τ))*(u⋅v)*ω^2)dΩ - ∫((K∘(τ))*divergence(u)*divergence(v))dΩ
      b(v) = ∫((v⋅nb)*p)dΓ

      op = AffineFEOperator(a, b, U, V)

      uh = solve(op)

      uh_y = (u->u[2])∘uh

      u = CellField(x->σ(x), Ω)

      
      # writevtk(Ω,mesh*".vtu", cellfields=["Real"=>real(uh_y),
      #                                                 "Imag"=>imag(uh_y),
      #                                                 "Real_Analitical"=>real(u),
      #                                                 "Imag_Analitical"=>imag(u),
      #                                                 "Error"=>abs(u-uh_y)])
      
      
      return 100 * sqrt(sum(∫(abs2(u - uh_y))*dΩ))/sqrt(sum(∫(abs2(u))*dΩ))
     
      

end


