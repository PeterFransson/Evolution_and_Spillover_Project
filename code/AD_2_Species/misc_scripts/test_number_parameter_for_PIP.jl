function create_syspar(R₀_aa_max::Real,
    γ::Real,
    σ²::Real,
    amplitude::Real,
    μ_a::Real,
    N_a::Integer,
    N_b::Integer,
    R₀_ratio::Real,
    Δᵣ::Real,
    c_ratio::Real)     
    #--Create system parameter-- 

    #maximum intraspecific basic reproduction number    
    R₀_bb_max = R₀_aa_max*R₀_ratio

    μ_b = μ_a+sqrt(2*σ²)*Δᵣ
    
    c_crit = min(R₀_aa_max,R₀_bb_max)*2/(R₀_aa_max+R₀_bb_max)
    c = c_crit*c_ratio      
            
    #Maximum intraspecific transmission rate
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
    β_ab_max = (β_aa_max+β_bb_max)/2*c

    syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude) 

    return syspar 
end 

function PIP_work_list()
    #Create syspar
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0 #variance and amplitude τ-function
    μ_a = 0.2 #Position of species a in resource space     
    N_a = 100000
    N_b = 10000  
    R₀_aa_max = 2.0 
    R₀_bb_max = 1.4  
    R₀_ratio = R₀_bb_max/R₀_aa_max #R₀_ratio∈(0,∞) 
    Δᵣ = 1.2 #Δᵣ∈(0,∞) Distance between species (1 Δᵣ = sqrt(2)*σ)
    c_ratio = 0.1 #c_ratio∈(0,1) IMPORTANT! This should not be to low otherwise its equal to a disconnected system 
    
    syspar = create_syspar(R₀_aa_max,
    γ,
    σ²,
    amplitude,
    μ_a,
    N_a,
    N_b,
    R₀_ratio,
    Δᵣ,
    c_ratio) 

    img_file_path = "./fig/misc/test_number_parameter_for_PIP/fig_4.svg"

    draw_PIP(img_file_path,syspar)    
end

PIP_work_list()