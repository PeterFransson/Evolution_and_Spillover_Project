function create_PIP(R₀_ratio::Real,Δᵣ::Real,c_ratio::Real,figure_name::AbstractString) 
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
    
    draw_PIP_n_coex_region(figure_name,syspar)
end 

function pip_work_list()
    R₀_ratio_vec = [1.0,1.5,2.0]
    Δᵣ_vec = collect(range(1.0,stop=3.0,length=4))
    c_ratio_vec = [0.8,0.6,0.2]

    for i in eachindex(R₀_ratio_vec)
        for j in eachindex(Δᵣ_vec)
            for k in eachindex(c_ratio_vec)
                R₀_ratio = R₀_ratio_vec[i]
                Δᵣ = Δᵣ_vec[i]
                c_ratio = c_ratio_vec[i]

                figure_name = "./fig/AD_2_Species/PIP_test/_R_$(i)_D_$(j)_C_$(k)_"

                create_PIP(R₀_ratio,Δᵣ,c_ratio,figure_name) 
            end            
        end        
    end
end

pip_work_list()