function σ(x)
    k_F = ω/c_F 
    k_P = ω/c_P
    τ = 1im * k_F
    γ = 1im * k_P
    α = ρ_P * c_P^2 * k_P / (ρ_F * c_F^2 * k_F)
    β = -P_0 / (ρ_F * c_F^2 * τ)
    ξ = exp(γ * t_P) + exp(-γ * t_P)
    ε = exp(-γ * t_P) - exp(γ * t_P)
    φ = exp((t_P - 2*(t_P+t_F)) * τ) + exp(-t_P * τ)

    # Porous Constants
    B_P = (β*(exp((t_P-(t_P+t_F))*τ) + 1/φ*(exp(-(t_P+t_F)*τ)-exp((2*t_P-3*(t_P+t_F))*τ)))) / (ε/φ*(exp(-t_P*τ) - exp((t_P-2*(t_P+t_F))*τ)) - α*ξ)
    A_P = -B_P
    # Fluid constants
    B_F = 1/φ * (B_P*ε - β*exp((t_P-(t_P+t_F))*τ))
    A_F = B_F*(exp(-2*(t_P+t_F)*τ)) + β*exp(-(t_P+t_F)*τ) 
    
    if (x[2] - t_P >= 0)
        return A_F * exp(1im * k_F * x[2]) + B_F * exp(-1im * k_F * x[2])
    else
        return  A_P * exp(1im * k_P * x[2]) + B_P * exp(-1im * k_P * x[2])
    end
end