# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Physical parameters
# Frequency parameters

function compute_ρP()
    return 1000
end

function compute_cP()
    return 4000
end

# Domains properties
ρ_F = 2000
c_F = 3000
ρ_P = compute_ρP()
c_P = compute_cP()

# Transducer pressure
P_0 = 5e5 + 5e5im
# Domain dimensions
L = 1 # Total length of the domain
t_P = 0.4 # Height of the porous domain
t_F = 1.6 # Height of the water domain
f = 25e3
ω = 2 * π * f