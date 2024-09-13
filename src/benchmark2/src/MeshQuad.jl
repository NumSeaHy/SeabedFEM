using Gmsh
using GridapGmsh
using Gridap.Io
using Gridap

include("Configuration.jl")

gmsh.initialize()
gmsh.model.add("fluid_and_porous_PML")



x(ρ, theta) =  ρ * cos(theta)
y(ρ, theta) = ρ * sin(theta)
DegToRad(x) = x * pi / 180.0
h = 0.1
p = 60
## POINTS

# Point Circles
cp = gmsh.model.geo.addPoint(0, 0, 0, h)
p1 = gmsh.model.geo.addPoint(r*x(r, DegToRad(0)), r*y(r, DegToRad(0)), 0, h)
p2 = gmsh.model.geo.addPoint(r*x(r, DegToRad(45)), r*y(r, DegToRad(45)), 0, h)
p3 = gmsh.model.geo.addPoint(r*x(r, DegToRad(90)), r*y(r, DegToRad(90)), 0, h)
p4 = gmsh.model.geo.addPoint(r*x(r, DegToRad(135)), r*y(r, DegToRad(135)), 0, h)
p5 = gmsh.model.geo.addPoint(r*x(r, DegToRad(180)), r*y(r, DegToRad(180)), 0, h)
p6 = gmsh.model.geo.addPoint(r*x(r, DegToRad(225)), r*y(r, DegToRad(225)), 0, h)
p7 = gmsh.model.geo.addPoint(r*x(r, DegToRad(270)), r*y(r, DegToRad(270)), 0, h)
p8 = gmsh.model.geo.addPoint(r*x(r, DegToRad(315)), r*y(r, DegToRad(315)), 0, h)
# Domain Points
p9 = gmsh.model.geo.addPoint(-L/2-d_PML, H/2 + d_PML , 0, h)
p10 = gmsh.model.geo.addPoint(-L/2, H/2 + d_PML, 0, h)
p11 = gmsh.model.geo.addPoint(0, H/2+d_PML, 0, h)
p12 = gmsh.model.geo.addPoint(L/2, H/2+d_PML, 0, h)
p13 = gmsh.model.geo.addPoint(L/2+d_PML, H/2+d_PML, 0, h)
p14 = gmsh.model.geo.addPoint(-L/2-d_PML, H/2, 0, h)
p15 = gmsh.model.geo.addPoint(-L/2, H/2, 0, h)
p16 = gmsh.model.geo.addPoint(0, H/2, 0, h)
p17 = gmsh.model.geo.addPoint(L/2, H/2, 0, h)
p18 = gmsh.model.geo.addPoint(L/2+d_PML, H/2, 0, h)
p19 = gmsh.model.geo.addPoint(-L/2-d_PML, 0, 0, h)
p20 = gmsh.model.geo.addPoint(-L/2, 0, 0, h)
p21 = gmsh.model.geo.addPoint(L/2, 0, 0, h)
p22 = gmsh.model.geo.addPoint(L/2+d_PML, 0, 0, h)
p23 = gmsh.model.geo.addPoint(-L/2-d_PML, -H/2, 0, h)
p24 = gmsh.model.geo.addPoint(-L/2, -H/2, 0, h)
p25 = gmsh.model.geo.addPoint(0, -H/2, 0, h)
p26 = gmsh.model.geo.addPoint(L/2, -H/2, 0, h)
p27 = gmsh.model.geo.addPoint(L/2+d_PML, -H/2, 0, h)
p28 = gmsh.model.geo.addPoint(-L/2-d_PML, -H/2-d_PML, 0, h)
p29 = gmsh.model.geo.addPoint(-L/2, -H/2-d_PML, 0, h)
p30 = gmsh.model.geo.addPoint(0, -H/2-d_PML, 0, h)
p31 = gmsh.model.geo.addPoint(L/2, -H/2-d_PML, 0, h)
p32 = gmsh.model.geo.addPoint(L/2+d_PML, -H/2-d_PML, 0, h)

# LINES
# Arcs
l1 = gmsh.model.geo.addCircleArc(p1, cp, p2)
l2 = gmsh.model.geo.addCircleArc(p2, cp, p3)
l3 = gmsh.model.geo.addCircleArc(p3, cp, p4)
l4 = gmsh.model.geo.addCircleArc(p4, cp, p5)
l5 = gmsh.model.geo.addCircleArc(p5, cp, p6)
l6 = gmsh.model.geo.addCircleArc(p6, cp, p7)
l7 = gmsh.model.geo.addCircleArc(p7, cp, p8)
l8 = gmsh.model.geo.addCircleArc(p8, cp, p1)
# Lines
l9 = gmsh.model.geo.addLine(p10, p9)
l10 = gmsh.model.geo.addLine(p15, p10) 
l11 = gmsh.model.geo.addLine(p14, p15)
l12 = gmsh.model.geo.addLine(p9, p14)
l13 = gmsh.model.geo.addLine(p11, p10)    
l14 = gmsh.model.geo.addLine(p15, p16)
l15 = gmsh.model.geo.addLine(p16, p11)
l16 = gmsh.model.geo.addLine(p12, p11)
l17 = gmsh.model.geo.addLine(p16, p17)
l18 = gmsh.model.geo.addLine(p17, p12)
l19 = gmsh.model.geo.addLine(p13, p12)
l20 = gmsh.model.geo.addLine(p17, p18)
l21 = gmsh.model.geo.addLine(p18, p13)
l22 = gmsh.model.geo.addLine(p14, p19)
l23 = gmsh.model.geo.addLine(p19, p20)
l24 = gmsh.model.geo.addLine(p20, p15)
l25 = gmsh.model.geo.addLine(p17, p21)
l26 = gmsh.model.geo.addLine(p21, p22)
l27 = gmsh.model.geo.addLine(p22, p18)
l28 = gmsh.model.geo.addLine(p19, p23)
l29 = gmsh.model.geo.addLine(p23, p24)
l30 = gmsh.model.geo.addLine(p24, p20)
l31 = gmsh.model.geo.addLine(p21, p26)
l32 = gmsh.model.geo.addLine(p26, p27)
l33 = gmsh.model.geo.addLine(p27, p22)
l34 = gmsh.model.geo.addLine(p23, p28)
l35 = gmsh.model.geo.addLine(p28, p29)
l36 = gmsh.model.geo.addLine(p29, p24)
l37 = gmsh.model.geo.addLine(p25, p24)
l38 = gmsh.model.geo.addLine(p29, p30)
l39 = gmsh.model.geo.addLine(p30, p25)
l40 = gmsh.model.geo.addLine(p26, p25)
l41 = gmsh.model.geo.addLine(p30, p31)
l42 = gmsh.model.geo.addLine(p31, p26)
l43 = gmsh.model.geo.addLine(p31, p32)
l44 = gmsh.model.geo.addLine(p32, p27)
l45 = gmsh.model.geo.addLine(p15, p4)
l46 = gmsh.model.geo.addLine(p16, p3)
l47 = gmsh.model.geo.addLine(p17, p2)
l48 = gmsh.model.geo.addLine(p21, p1)
l49 = gmsh.model.geo.addLine(p26, p8)
l50 = gmsh.model.geo.addLine(p25, p7)
l51 = gmsh.model.geo.addLine(p24, p6)
l52 = gmsh.model.geo.addLine(p20, p5)

source = gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5,l6, l7, l8])
gmsh.model.setPhysicalName(1, source, "source")
l_pml = gmsh.model.addPhysicalGroup(1, [l9, l12, l22, l28, l34, l35, l38, l41, l43, l44, l33, l27, l21, l19, l16, l13])
gmsh.model.setPhysicalName(1, l_pml, "sides")
# CURVE lOOPS
# Interior sectors
cl1 = gmsh.model.geo.addCurveLoop([-l25, l47, -l1, -l48])
cl2 = gmsh.model.geo.addCurveLoop([-l47, -l17, l46, -l2])
cl3 = gmsh.model.geo.addCurveLoop([-l46, -l14, l45, -l3])
cl4 = gmsh.model.geo.addCurveLoop([-l45, -l24, l52, -l4])
cl5 = gmsh.model.geo.addCurveLoop([-l52, -l30, l51, -l5])
cl6 = gmsh.model.geo.addCurveLoop([-l51, -l37, l50, -l6])
cl7 = gmsh.model.geo.addCurveLoop([-l50, -l40, l49, -l7])
cl8 = gmsh.model.geo.addCurveLoop([-l49, -l31, l48, -l8])
# Exterior
cl9 = gmsh.model.geo.addCurveLoop([l9, l12, l11, l10])
cl10 = gmsh.model.geo.addCurveLoop([l13, -l10, l14, l15])
cl11 = gmsh.model.geo.addCurveLoop([l16, -l15, l17, l18])
cl12 = gmsh.model.geo.addCurveLoop([l19, -l18, l20, l21])
cl13 = gmsh.model.geo.addCurveLoop([-l11, l22, l23, l24])
cl14 = gmsh.model.geo.addCurveLoop([-l20, l25, l26, l27])
cl15 = gmsh.model.geo.addCurveLoop([-l23, l28, l29, l30])
cl16 = gmsh.model.geo.addCurveLoop([-l26, l31, l32, l33])
cl17 = gmsh.model.geo.addCurveLoop([-l29, l34, l35, l36])
cl18 = gmsh.model.geo.addCurveLoop([l37, -l36, l38, l39])
cl19 = gmsh.model.geo.addCurveLoop([l40, -l39, l41, l42])
cl20 = gmsh.model.geo.addCurveLoop([-l32, -l42, l43, l44])


lines = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10,
        l11, l12, l13, l14, l15, l16, l17, l18, l19, l20,
        l21, l22, l23, l24, l25, l26, l27, l28, l29, l30,
        l31, l32, l33, l34, l35, l36, l37, l38, l39, l40,
        l41, l42, l43, l44, l45, l46, l47, l48, l49, l50,
        l51, l52]

for i in lines
    gmsh.model.geo.mesh.setTransfiniteCurve(i, p)
end

# SURFACES
# Interior sectors
ps1 = gmsh.model.geo.addPlaneSurface([cl1])
ps2 = gmsh.model.geo.addPlaneSurface([cl2])
ps3 = gmsh.model.geo.addPlaneSurface([cl3])
ps4 = gmsh.model.geo.addPlaneSurface([cl4])
ps5 = gmsh.model.geo.addPlaneSurface([cl5])
ps6 = gmsh.model.geo.addPlaneSurface([cl6])
ps7 = gmsh.model.geo.addPlaneSurface([cl7])
ps8 = gmsh.model.geo.addPlaneSurface([cl8])
# Exterior
ps9 = gmsh.model.geo.addPlaneSurface([cl9])
ps10 = gmsh.model.geo.addPlaneSurface([cl10])
ps11 = gmsh.model.geo.addPlaneSurface([cl11])
ps12 = gmsh.model.geo.addPlaneSurface([cl12])
ps13 = gmsh.model.geo.addPlaneSurface([cl13])
ps14 = gmsh.model.geo.addPlaneSurface([cl14])
ps15 = gmsh.model.geo.addPlaneSurface([cl15])
ps16 = gmsh.model.geo.addPlaneSurface([cl16])
ps17 = gmsh.model.geo.addPlaneSurface([cl17])
ps18 = gmsh.model.geo.addPlaneSurface([cl18])
ps19 = gmsh.model.geo.addPlaneSurface([cl19])
ps20 = gmsh.model.geo.addPlaneSurface([cl20])

surfaces = gmsh.model.addPhysicalGroup(2, [ps1, ps2, ps3, ps4, ps5, ps6, ps7, ps8, ps9, ps10, ps11, ps12, ps13, ps14, ps15, ps16, ps17, ps18, ps19, ps20])
gmsh.model.setPhysicalName(2, surfaces, "physicalDomain")

# for i in lines
#     gmsh.model.geo.mesh.setTransfiniteCurve(i, 60)
# end
surfaces = [ps1, ps2, ps3, ps4, ps5, ps6, ps7, ps8, ps9, ps10, ps11, ps12, ps13, ps14, ps15, ps16, ps17, ps18, ps19, ps20]
for i in surfaces
    gmsh.model.geo.mesh.setTransfiniteSurface(i)
    gmsh.model.geo.mesh.setRecombine(2, i)
end



# Sincronización y generación de la malla
gmsh.model.geo.synchronize()


gmsh.model.mesh.generate(2)

# Guardado de la malla y finalización
gmsh.write("./data/Qudrangular.msh")
gmsh.finalize()

# model = GmshDiscreteModel("./data/coarse2.msh")

# fn = "./data/pml_mesh_quad.json"
# to_json_file(model,fn)

# writevtk(model, "results/coarse2")