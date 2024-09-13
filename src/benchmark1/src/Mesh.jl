using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

# Input the physical parameters of the problem to calculate meshes sizes
include("Configuration.jl")

gmsh.initialize()
gmsh.model.add("fluid_and_porous")

# Determining the mesh size
h1 = c_F/(15*f)
h2 = c_P/(15*f)

# Corrners of the 2 Domains
gmsh.model.geo.addPoint(0, 0, 0, h2, 1)
gmsh.model.geo.addPoint(L, 0, 0, h2, 2)
gmsh.model.geo.addPoint(L, t_F+t_P, 0, h1, 3)
gmsh.model.geo.addPoint(0, t_F+t_P, 0, h1, 4)
gmsh.model.geo.addPoint(L, t_P, 0, min(h1, h2), 5)
gmsh.model.geo.addPoint(0, t_P, 0, min(h1, h2), 6)

# Edges of the 2 Domains
gmsh.model.geo.add_line(1, 2, 1)
gmsh.model.geo.add_line(2, 5, 2)
gmsh.model.geo.add_line(5, 3, 3)
gmsh.model.geo.add_line(3, 4, 4)
gmsh.model.geo.add_line(4, 6, 5)
gmsh.model.geo.add_line(6, 1, 6)
gmsh.model.geo.add_line(5, 6, 7) # Interface Line

# Curve loop and surface
gmsh.model.geo.add_curve_loop([1, 2, 7, 6], 1)
gmsh.model.geo.add_curve_loop([3, 4, 5, -7], 2)

gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([2], 2)

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
gmsh.model.addPhysicalGroup(1, [1], 1)
gmsh.model.addPhysicalGroup(1, [2], 2)
gmsh.model.addPhysicalGroup(1, [3], 3)
gmsh.model.addPhysicalGroup(1, [4], 4)
gmsh.model.addPhysicalGroup(1, [5], 5)
gmsh.model.addPhysicalGroup(1, [6], 6)
gmsh.model.addPhysicalGroup(1, [7], 7)
# Set physical groups for 2D entities
gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.addPhysicalGroup(2, [2], 2)
# Set physical names for 1D and 2D entities
gmsh.model.setPhysicalName(1, 1, "bottom")
gmsh.model.setPhysicalName(1, 2, "right_porous")
gmsh.model.setPhysicalName(1, 3, "right_water")
gmsh.model.setPhysicalName(1, 4, "transducer")
gmsh.model.setPhysicalName(1, 4, "transducer")
gmsh.model.setPhysicalName(1, 5, "left_water")
gmsh.model.setPhysicalName(1, 6, "left_porous")
gmsh.model.setPhysicalName(1, 7, "interface")
gmsh.model.setPhysicalName(2, 1, "porous_domain")
gmsh.model.setPhysicalName(2, 2, "water_domain")

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/benchmark_2_coarse.msh")

gmsh.finalize()

# This steps are to visualize the mesh and his labels using Gridap's API and to write the 
# mesh to .json extension. Extension used by Gridap to run the cases.
model = GmshDiscreteModel("./data/benchmark_2_coarse.msh")

writevtk(model,"./data/benchmark_2_coarse")

fn = "./data/benchmark_2_coarse.json"
to_json_file(model,fn)