struct parameters
    R₀_aa::Real #Maximum Basic reproduction number in species a
    R₀_bb::Real #Maximum Basic reproduction number in species b
    R₀_ab::Real #Maximum intraspecific basic reproduction number
    N_a::Integer #Pupulation size of species a
    N_b::Integer #Pupulation size of species b
end

function create_syspar(R₀_ratio::Real,Δᵣ::Real,c_ratio::Real)     
    #--Create system parameter-- 

    #maximum intraspecific basic reproduction number
    R₀_aa_max = 1.4
    R₀_bb_max = R₀_aa_max*R₀_ratio

    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2
    μ_b = μ_a+sqrt(2*σ²)*Δᵣ
    
    c_crit = min(R₀_aa_max,R₀_bb_max)*2/(R₀_aa_max+R₀_bb_max)
    c = c_crit*c_ratio
        
    N_a = 10^3
    N_b = 10^3    
            
    #Maximum intraspecific transmission rate
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
    β_ab_max = (β_aa_max+β_bb_max)/2*c

    syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude) 

    return syspar 
end 
