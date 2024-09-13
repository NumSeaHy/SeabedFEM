"""
This file generates the mesh of the PML domain using Gmsh 
The mesh includes the fluid domain, the porous domain and the PML and it is centered at the origin.
The mesh is saved in the folder data in the format .msh and .json
"""
# Load the packages
using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

# Input of the physical parameters neccesary to calculate the mesh size
include("Configuration.jl")

# Initialize gmsh
gmsh.initialize()
gmsh.model.add("fluid_and_porous_PML")

# Mesh size depending on the angular frequency
h1 = 2 * π * real(c_F(ω))/(15*ω) 
h2 = 2 * π * real(c_P(ω))/(15*ω)
h_clam = h1/6


# Function to use polar coordinates
x(r, theta) =  r * cos(theta)
y(r, theta) = r * sin(theta)

# Corners of Porous Domain
p_P1 = gmsh.model.geo.addPoint(-L/2, 0, 0, h2)
p_P2 = gmsh.model.geo.addPoint(L/2, 0, 0, h2)
p_P3 = gmsh.model.geo.addPoint(L/2, t_P, 0, min(h1, h2)) 
p_P4 = gmsh.model.geo.addPoint(-L/2, t_P, 0, min(h1, h2))

# Corners of the Porous PML
    # Left Part
    p_PML_P1 = gmsh.model.geo.addPoint(-d_PML-L/2, 0, 0, h2)
    p_PML_P2 = gmsh.model.geo.addPoint(-d_PML-L/2, t_P, 0, min(h1, h2)) # == p_PML_F1
    # Right Part
    p_PML_P3 = gmsh.model.geo.addPoint(L/2+d_PML, 0, 0, h2)
    p_PML_P4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P, 0, min(h1, h2)) # == p_PML_F3 

# Corrners of Fluid Domain
p_F1 = gmsh.model.geo.addPoint(L/2, t_P+t_F, 0, h1) 
p_F2 = gmsh.model.geo.addPoint(-L/2, t_P+t_F, 0, h1)

# Corrners of the Transducer
p_T1 = gmsh.model.geo.addPoint(t/2, t_P+t_F, 0, h1)
p_T2 = gmsh.model.geo.addPoint(-t/2, t_P+t_F, 0, h1)

# Corners of the Fluid PML
    # Left Part
    #p_PML_F1 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P, 0, h1)
    p_PML_F2 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P+t_F, 0, h1)
    # Right Part
    #p_PML_F3 = gmsh.model.geo.addPoint(L/2+d_PML, t_P, 0, h1)
    p_PML_F4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P+t_F, 0, h1)

# Corners of the Bottom PML
p_PML_B1 = gmsh.model.geo.addPoint(-L/2, -d_PML, 0, h2)
p_PML_B2 = gmsh.model.geo.addPoint(-L/2-d_PML, -d_PML, 0, h2)
#p_PML_B3 = gmsh.model.geo.addPoint(-L/2-d_PML, 0, 0, h2)
p_PML_B4 = gmsh.model.geo.addPoint(L/2, -d_PML, 0, h2)
p_PML_B5 = gmsh.model.geo.addPoint(L/2+d_PML, -d_PML, 0, h2)
#p_PML_B6 = gmsh.model.geo.addPoint(L/2+d_PML, 0, 0, h2)

# Clam points
centerPoint = gmsh.model.geo.addPoint(x_c, y_c, 0, h_clam)
p2 = gmsh.model.geo.addPoint(x_c + x(r_in, θ_b + θ_o/2), y_c + y(r_in, θ_b + θ_o/2), 0, h_clam)
p3 = gmsh.model.geo.addPoint(x_c + x(r_in, θ_b + θ_o/2 + pi), y_c + y(r_in, θ_b + θ_o/2 + pi), 0, h_clam)
p4 = gmsh.model.geo.addPoint(x_c + x(r_in, θ_b - θ_o/2), y_c + y(r_in, θ_b - θ_o/2), 0, h_clam)
p5 = gmsh.model.geo.addPoint(x_c + x(r_out, θ_b + θ_o/2), y_c + y(r_out, θ_b + θ_o/2), 0, h_clam)
p6 = gmsh.model.geo.addPoint(x_c + x(r_out, θ_b + θ_o/2 + pi), y_c + y(r_out, θ_b + θ_o/2 + pi), 0, h_clam)
p7 = gmsh.model.geo.addPoint(x_c + x(r_out, θ_b - θ_o/2), y_c + y(r_out, θ_b - θ_o/2), 0, h_clam)

# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p_P1, p_P2)
l_P2 = gmsh.model.geo.add_line(p_P2, p_P3)
l_P3 = gmsh.model.geo.add_line(p_P3, p_P4)
l_P4 = gmsh.model.geo.add_line(p_P4, p_P1)

# Edges of the Porous PML
    # Left Part
    l_PML_P1 = -l_P4
    l_PML_P2 = gmsh.model.geo.add_line(p_P4, p_PML_P2)
    l_PML_P3 = gmsh.model.geo.add_line(p_PML_P2, p_PML_P1)
    l_PML_P4 = gmsh.model.geo.add_line(p_PML_P1, p_P1)
    # Right Part
    l_PML_P5 = -l_P2
    l_PML_P6 = gmsh.model.geo.add_line(p_P2, p_PML_P3)
    l_PML_P7 = gmsh.model.geo.add_line(p_PML_P3, p_PML_P4)
    l_PML_P8 = gmsh.model.geo.add_line(p_PML_P4, p_P3)

# Edges of the Fluid
l_F1 = -l_P3
l_F2 = gmsh.model.geo.add_line(p_P3, p_F1)
l_T1 = gmsh.model.geo.add_line(p_F1, p_T1)
l_T2 = gmsh.model.geo.add_line(p_T1, p_T2)
l_F3 = gmsh.model.geo.add_line(p_T2, p_F2)
l_F4 = gmsh.model.geo.add_line(p_F2, p_P4)

# Edges of the Fluid PML
    # Left Part
    l_PML_F1 = -l_F4
    l_PML_F2 = gmsh.model.geo.add_line(p_F2, p_PML_F2)
    l_PML_F3 = gmsh.model.geo.add_line(p_PML_F2, p_PML_P2)
    l_PML_F4 = -l_PML_P2

    # Right Part
    l_PML_F5 = -l_F2
    l_PML_F6 = -l_PML_P8
    l_PML_F7 = gmsh.model.geo.add_line(p_PML_P4, p_PML_F4)
    l_PML_F8 = gmsh.model.geo.add_line(p_PML_F4, p_F1)

# Edges of the Bottom PML
    # Left Part
    l_PML_B1 = gmsh.model.geo.add_line(p_PML_B2, p_PML_B1)
    l_PML_B2 = gmsh.model.geo.add_line(p_PML_B1, p_P1)
    l_PML_B3 = -l_PML_P4
    l_PML_B4 = gmsh.model.geo.add_line(p_PML_P1, p_PML_B2)
    # Center Part
    l_PML_B5 = gmsh.model.geo.add_line(p_PML_B1, p_PML_B4)
    l_PML_B6 = gmsh.model.geo.add_line(p_PML_B4, p_P2)
    l_PML_B7 = -l_P1
    l_PML_B8 = -l_PML_B2 # gmsh.model.geo.add_line(p_P1, p_PML_B1)
    # Right Part
    l_PML_B9 = gmsh.model.geo.add_line(p_PML_B4, p_PML_B5)
    l_PML_B10 = gmsh.model.geo.add_line(p_PML_B5, p_PML_P3)
    l_PML_B11 = -l_PML_P6
    l_PML_B12 = -l_PML_B6 # gmsh.model.geo.add_line(p_P2, p_PML_B4)

# Edges of the clam
arc_i1 = gmsh.model.geo.addCircleArc(p2, centerPoint, p3)
arc_i2 = gmsh.model.geo.addCircleArc(p3, centerPoint, p4)
line_1 = gmsh.model.geo.addLine(p4, p7)
arc_e1 = gmsh.model.geo.addCircleArc(p6, centerPoint, p7)
arc_e2 = gmsh.model.geo.addCircleArc(p5, centerPoint, p6)
line_2 = gmsh.model.geo.addLine(p5, p2)

# Curve Loops
cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])
cl_PML_P1 = gmsh.model.geo.add_curve_loop([l_PML_P1, l_PML_P2, l_PML_P3, l_PML_P4])
cl_PML_P2 = gmsh.model.geo.add_curve_loop([l_PML_P5, l_PML_P6, l_PML_P7, l_PML_P8])
cl_F = gmsh.model.geo.add_curve_loop([l_F1, l_F2, l_T1, l_T2, l_F3, l_F4])
cl_PML_F1 = gmsh.model.geo.add_curve_loop([l_PML_F1, l_PML_F2, l_PML_F3, l_PML_F4])
cl_PML_F2 = gmsh.model.geo.add_curve_loop([l_PML_F5, l_PML_F6, l_PML_F7, l_PML_F8])
cl_PML_B1 = gmsh.model.geo.add_curve_loop([l_PML_B1, l_PML_B2, l_PML_B3, l_PML_B4])
cl_PML_B2 = gmsh.model.geo.add_curve_loop([l_PML_B5, l_PML_B6, l_PML_B7, l_PML_B8])
cl_PML_B3 = gmsh.model.geo.add_curve_loop([l_PML_B9, l_PML_B10, l_PML_B11, l_PML_B12])
cl_rock = gmsh.model.geo.addCurveLoop([arc_i1, arc_i2, line_1, -arc_e1, -arc_e2, line_2])

# Surfaces
s_P = gmsh.model.geo.addPlaneSurface([cl_P, cl_rock])
s_PML_P1 = gmsh.model.geo.addPlaneSurface([cl_PML_P1])
s_PML_P2 = gmsh.model.geo.addPlaneSurface([cl_PML_P2])
s_F = gmsh.model.geo.addPlaneSurface([cl_F])
s_PML_F1 = gmsh.model.geo.addPlaneSurface([cl_PML_F1])
s_PML_F2 = gmsh.model.geo.addPlaneSurface([cl_PML_F2])
s_PML_B1 = gmsh.model.geo.addPlaneSurface([cl_PML_B1])
s_PML_B2 = gmsh.model.geo.addPlaneSurface([cl_PML_B2])
s_PML_B3 = gmsh.model.geo.addPlaneSurface([cl_PML_B3])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
f_T = gmsh.model.addPhysicalGroup(1, [l_PML_F8, l_T1, l_T2, l_F3, l_PML_F2])
f_N = gmsh.model.addPhysicalGroup(1, [ l_PML_F3, l_PML_P3, l_PML_B4, l_PML_B10, l_PML_P7,  l_PML_F7])
f_B = gmsh.model.addPhysicalGroup(1, [l_PML_B1, l_PML_B5, l_PML_B9])
l_rock = gmsh.model.addPhysicalGroup(1, [arc_i1, arc_i2, line_1, arc_e1, arc_e2, line_2])

# Set physical groups for 2D entities
# f_P = gmsh.model.addPhysicalGroup(2, [s_P])
f_PML_P = gmsh.model.addPhysicalGroup(2, [s_PML_P1, s_PML_P2])
# f_F = gmsh.model.addPhysicalGroup(2, [s_F])
f_PML_F = gmsh.model.addPhysicalGroup(2, [s_PML_F1, s_PML_F2])
f_PML_B_s = gmsh.model.addPhysicalGroup(2, [s_PML_B1, s_PML_B3])
f_PML_B = gmsh.model.addPhysicalGroup(2, [s_PML_B2])
f_Ph = gmsh.model.addPhysicalGroup(2, [s_F, s_P])
# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, f_T, "transducer")
gmsh.model.setPhysicalName(1, f_N, "sides")
gmsh.model.setPhysicalName(1, f_B, "bottom")
gmsh.model.setPhysicalName(1, l_rock, "rock")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_Ph, "physical_domain")

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/mesh.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/mesh.msh")

# Write the mesh to a vtk file
writevtk(model,"./results/mesh")
