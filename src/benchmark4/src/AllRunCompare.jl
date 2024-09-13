"""
This script calls sequantially all the scripts to run the simulations and compare the solutions between the conventional approach and the translated approach.
"""

include("./Mesh.jl")
include("./Run.jl")
include("./RunTranslated.jl")
include("./PostProcess.jl")
