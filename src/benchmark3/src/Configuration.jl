# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 1
t_P = 0.5 # Heigth of the porous domain
t_F = 1.5 # Heigth of the fluid domain
t = 0.25 # Length of the sonar
# d_PML = 0.125
d_PML = 0.1

# Check possible Domain
if t > L
    error("The length of the sonar cannot be bigger than the length of the domain")
end

# Frequency parameters
f = 30e3
ω = 2 * π * f

function compute_ρP(ω)
    return 1000
end

function compute_cP(ω)
    return 10e3
end

# Domains properties
ρ_F = 2000
c_F = 8e3
ρ_P = 1000
c_P = 6e3

# Transducer pressure
P_0 = 5e5 + 5e5im

RPML = 10
σ_0 = -3/4 * log(RPML)/d_PML
