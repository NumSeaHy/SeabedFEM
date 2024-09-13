using Gridap
using Gridap.Geometry
using Gridap.Fields
using Gridap.Io
using GridapGmsh
using JLD2

include("./FunMesh.jl")
include("./Configuration.jl")

function Run(f, save_vtk::Bool, save_jld::Bool)
    
    L2 = zeros(length(f))
    
    for (i, data) in enumerate(f)
      
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
    
        name = "data_freq/mesh"*string(data/1e3)*"KHz"
         
        MeshGenerator(c_F, c_P, data, L, t_P+t_F, t_P, name)
        
        model = DiscreteModelFromFile(name*".json")
    
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
        ω = 2 * π * data
        p = P_0
        a(u, v) = ∫((ρ∘(τ))*(u⋅v)*ω^2)dΩ - ∫((K∘(τ))*divergence(u)*divergence(v))dΩ
        b(v) = ∫((v⋅nb)*p)dΓ
    
        op = AffineFEOperator(a, b, U, V)
    
        uh = solve(op)
        
        u_x = (u->u[1])∘uh
        u_y = (u->u[2])∘uh

        if save_vtk
            save_name = "results_freq/peak$(i)"
            writevtk(Ω,save_name*".vtu", cellfields=["Re_x"=>real(u_x),
                                                     "Im_x"=>imag(u_x),
                                                     "Re_y"=>real(u_y),                                            
                                                     "Im_y"=>imag(u_y)])
        end 
    
        L2[i] = sqrt(sum(∫(abs2((u_x)))*dΩ) + sum(∫(abs2((u_y)))*dΩ))         
    end
    if save_jld
        @save "results_freq/freq.jld2" f L2
    end     
end

