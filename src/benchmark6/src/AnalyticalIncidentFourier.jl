module AnalyticalIncidentFourier

using Gridap
using SpecialFunctions
using FFTW


export  EndFireParams, u_incident_modes, compute_scat_coefs

struct EndFireParams
    L::Float64
    t_P::Float64
    t_F::Float64
    d_PML::Float64
    Nₛ::Int
    xᵦ::Float64
    yᵦ::Float64
    Δy::Float64
    k_F::Float64
    k_P::ComplexF64
    σ_0::Float64
    pml_active::Bool
    ρ_F::Float64
    ρ_P:: ComplexF64
    ω::Float64
    P0::Float64  
end

"""
Strecthing function for the PML layer on the interface.
    
"""
x̂(x, L, k_F, σ_0, d_PML, pml_active) = pml_active ? (x >= L/2 ? x + 1im/k_F*σ_0*1/(d_PML^2) * (x-L/2)^3 : x <= -L/2 ? x + 1im/k_F*σ_0*1/(d_PML^2) * (x+L/2)^3 : x) : x

function π_end_fire_array(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i*params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    πF = 1im/4 * sum(A[j] * hankelh1(0, params.k_F * sqrt((x̂(x, params.L, params.k_F, params.σ_0, params.d_PML, params.pml_active) .- xpos[j])^2 + (params.t_P - ypos[j])^2)) for j in eachindex(ypos))
    
    return πF
end

function u_end_fire_array(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i * params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    uF = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x̂(x, params.L, params.k_F, params.σ_0, params.d_PML, params.pml_active) .- xpos[j])^2 + (params.t_P - ypos[j])^2)) * (params.t_P - ypos[j]) / sqrt((x̂(x, params.L, params.k_F, params.σ_0, params.d_PML, params.pml_active) .- xpos[j])^2 + (params.t_P - ypos[j])^2) for j in eachindex(ypos))
    
    return uF
end

function compute_fft_coefs(N, func, params::EndFireParams)
    
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / N # Grid step size
    x_points = range(-a, stop=a - h_fft, length=N) # Discretization of the interface domain [m], with size N

    # Compute the wavenumbers
    k = 2 * π * fftfreq(N, 1) ./ h_fft

    # Compute the values of the function on the grid points
    f_values = map(x -> func(x, params), x_points)

    # Compute the Fourier coefficients of the input data and normalize by N using the FFT
    coefs = fft(f_values) ./ N

    return coefs, k
end

function compute_scat_coefs(N, params::EndFireParams)
    
    π0, k = compute_fft_coefs(N, π_end_fire_array, params)
    u0, _ = compute_fft_coefs(N, u_end_fire_array, params)

    
    π0[1] = params.P0 # Incident pressure field
    u0[1] = -1im * params.k_F * params.P0/(params.ρ_F*params.ω^2) # Incident displacement field

    Βsel_F(k_F, k) = k_F > abs(k) ? 1im.*sqrt(k_F^2 - k^2) : -sqrt(k^2 - k_F^2)
    Βsel_P(k_P, k) = abs(k_P) > abs(k) ? -1im.*sqrt(k_P^2 - k^2) : +sqrt(k^2 - k_P^2)
    
    βF = Βsel_F.(params.k_F, k) # Bethas in the fluid domain
    βP = Βsel_P.(params.k_P, k) # Bethas in the porous domain
    
    # Coefficients of the 2x2 linear system associated with the coupling conditions on the interface
    A = 1 / (params.ρ_F * params.ω^2) .* βF
    B = 1 / (params.ρ_P * params.ω^2) .* βP # I have changed the sign of the B coefficient

    # Solution of the 2x2 linear system
    πF_s = (B .* π0 .- u0) ./ (A .- B)
    πP = πF_s .+ π0

    return πF_s, πP, βF, βP, k
end

function uF_s(x, a, t_P, βF, πF_s, k, ρ_F, ω)
    
    ux = 1 / (ρ_F * ω^2) * sum(πF_s[j] * 1im * k[j] * exp(βF[j] * (x[2] - t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βF))
    uy = 1 / (ρ_F * ω^2) * sum(πF_s[j] * βF[j] * exp(βF[j] * (x[2] - t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βF))
    
    return VectorValue(ux, uy)
end

function uP(x, a, t_P, βP, πP, k, ρ_P, ω)
    
    ux = 1 / (ρ_P * ω^2) * sum(πP[j] * 1im * k[j] * exp(βP[j] * (x[2] - t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βP))
    uy = 1 / (ρ_P * ω^2) * sum(πP[j] * βP[j] * exp(βP[j] * (x[2] - t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βP))
    
    return VectorValue(ux, uy)
end

function u_incident_modes(x, πF_s, πP, βF, βP, k,  params::EndFireParams)    
    
    if abs(x[1]) <= params.L/2 && (x[2] - params.t_P) < 0 && x[2] >= 0
        return uP(x, params.L/2 + params.d_PML, params.t_P, βP, πP, k, params.ρ_P, params.ω) 
    elseif abs(x[1]) <= params.L/2 && (x[2] - params.t_P) >= 0 && (x[2] - (params.t_P + params.t_F)) < 0
        return uF_s(x, params.L/2 + params.d_PML, params.t_P, βF, πF_s, k, params.ρ_F, params.ω)
    else
        return VectorValue(0.0, 0.0)
    end
end

"""
This function implements the analytical solution for the incident field in the two domain problem using the Sommerfeld boundary condition.
"""
function exact_solution_som(x, ρ_F, c_F, k_F, ρ_P, c_P, k_P, P_0, t_F, t_P, d_PML, σ_0)

    # Characteristic impedance at the fluid and porous domains
    Z_F = ρ_F * c_F
    Z_P = ρ_P * c_P

    # Source term
    α = P_0/(1im*ω*Z_F)
    
    # Fluid amplitudes
    A_F = (Z_F-Z_P)/Z_F * α * exp(1im*t_F*k_F)/((Z_F+Z_P)/Z_F*exp(1im*t_P*k_F)+(Z_P-Z_F)/Z_F*exp(1im*(2*t_F+t_P)*k_F))
    B_F = A_F*exp(2im*(t_F+t_P)*k_F) + α*exp(1im*(t_F+t_P)*k_F)
    
    # Porous amplitudes
    B_P = (A_F*(exp(1im*t_P*k_F)+exp(1im*(2*t_F+t_P)*k_F))+α*exp(1im*t_F*k_F))*exp(1im*t_P*k_P) 
   
    # Evaluation at fluid, porous and bottom PML
    if (abs(x[1]) <= L/2 && x[2] - t_P >= 0 && (x[2]-(t_P+t_F)) < 0)
        return VectorValue(0., A_F * exp(1im * k_F * x[2]) + B_F * exp(-1im * k_F * x[2])) # fluid domain
    elseif (abs(x[1]) <= L/2 && x[2] - t_P < 0 && x[2]>0)
        return VectorValue(0., B_P * exp(-1im * k_P * x[2])) # porous domain
    else
        # return VectorValue(0., B_P * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2)) # bottom PML (x[2]<0)
        return VectorValue(0., 0.) # bottom PML (x[2]<0)
    end
end

function sommerfeld_scat(x, ρ_F, c_F, k_F, ρ_P, c_P, k_P, P_0, t_F, t_P)

    # Characteristic impedance at the fluid and porous domains
    Z_F = ρ_F * c_F
    Z_P = ρ_P * c_P

    # Source term
    α = P_0/(1im*ω*Z_F)
    
    # Fluid amplitudes
    A_F = (Z_F-Z_P)/Z_F * α * exp(1im*t_F*k_F)/((Z_F+Z_P)/Z_F*exp(1im*t_P*k_F)+(Z_P-Z_F)/Z_F*exp(1im*(2*t_F+t_P)*k_F))
    B_F = A_F*exp(2im*(t_F+t_P)*k_F) + α*exp(1im*(t_F+t_P)*k_F)
    
    # Porous amplitudes
    B_P = (A_F*(exp(1im*t_P*k_F)+exp(1im*(2*t_F+t_P)*k_F))+α*exp(1im*t_F*k_F))*exp(1im*t_P*k_P) 
   
    # Evaluation at fluid, porous and bottom PML
    if (abs(x[1]) <= L/2 && x[2] - t_P >= 0 && (x[2]-(t_P+t_F)) < 0)
        return VectorValue(0., A_F * exp(1im * k_F * x[2])) # fluid domain
    elseif (abs(x[1]) <= L/2 && x[2] - t_P < 0 && x[2]>0)
        return VectorValue(0., B_P * exp(-1im * k_P * x[2])) # porous domain
    else
        # return VectorValue(0., B_P * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2)) # bottom PML (x[2]<0)
        return VectorValue(0., 0.) # bottom PML (x[2]<0)
    end
end


"""
This function implements the calculation of the coefficients of the analytical solution for the incident field in the two domain problem using the
Perfectly Matched Layer (PML) boundary condition.
"""
    function solve_coefficients_pml(ω, Z_F, Z_P, P_0, t_P, t_F, k_F, k_P, σ_0, d_PML)
        
        # Matrix of the system
        M = Complex{Float64}[
            exp(1im*(t_P+t_F)*k_F) -exp(-1im*(t_P+t_F)*k_F) 0 0 0 0;
            exp(1im*t_P*k_F) exp(-1im*t_P*k_F) -exp(1im*t_P*k_P) -exp(-1im*t_P*k_P) 0 0;
            exp(1im*t_P*k_F) -exp(-1im*t_P*k_F) -Z_P/Z_F*exp(1im*t_P*k_P) Z_P/Z_F*exp(-1im*t_P*k_P) 0 0;
            0 0 1 1 -1 -1; 
            0 0 1 -1 -1 1;
            # 0 0 0 0 1 0]  
            0 0 0 0 exp(-1im * k_P * d_PML)*exp(d_PML*σ_0/3) exp(1im * k_P * d_PML)*exp(-d_PML*σ_0/3)]

        # Independent term of the linear system
        b = Complex{Float64}[
            -P_0 / (1im*ω*Z_F),
            0,
            0,
            0,  
            0,  
            0   
        ]

        
        solution = M \ b

        A_F = solution[1]
        B_F = solution[2]
        A_P = solution[3]
        B_P = solution[4]
        A_PML = solution[5]
        B_PML = solution[6]

        return A_F, B_F, A_P, B_P, A_PML, B_PML
    end

    function exact_quadratic_pml(x, t_P, k_F, k_P, d_PML, σ_0, A_F, B_F, A_P, B_P, A_PML, B_PML)
        if (abs(x[1]) <= L/2 && x[2] - t_P >= 0 && (x[2]-(t_P+t_F)) < 0)
            return VectorValue(0., A_F * exp(1im * k_F * x[2]) + B_F * exp(-1im * k_F * x[2]))
        elseif (abs(x[1]) <= L/2 && x[2] - t_P < 0 && x[2]>0)
            return  VectorValue(0., A_P * exp(1im * k_P * x[2]) + B_P * exp(-1im * k_P * x[2]))
        else
            # return  VectorValue(0., A_PML * exp(1im * k_P * x[2]) * exp(-σ_0/3*(x[2]^3)/d_PML^2) + B_PML * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2))
            return  VectorValue(0., 0.)
        end
    end

end 