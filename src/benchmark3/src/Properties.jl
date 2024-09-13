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