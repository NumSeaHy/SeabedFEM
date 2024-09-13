using Gmsh
using Gridap

"""
    refine_mesh(input_mesh, output_mesh)

This function takes an input mesh generated with Gmsh and refines it, producing a refined version of the 
input mesh. This generates an interlocked mesh that preserves the original nodes, allowing for an 
efficient mesh convergence study.

# Arguments
- `input_mesh`: The mesh to be refined.
- `output_mesh`: The resulting refined mesh.

# Example
```julia
refine_mesh("path/to/mesh1", "path/to/mesh2")
"""
function refine_mesh(input_mesh, output_mesh)
    gmsh.initialize()
    gmsh.open(input_mesh*".msh")

    gmsh.model.mesh.refine()

    gmsh.write(output_mesh*".msh")
    
    gmsh.finalize()
    
    # model = GmshDiscreteModel(output_mesh*".msh")
    # writevtk(model, output_mesh)
    # to_json_file(model, output_mesh*".json")
    
end

