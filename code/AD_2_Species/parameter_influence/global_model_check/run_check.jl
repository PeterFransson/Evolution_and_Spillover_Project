function compute() 
    #we assume recovery rate (γ) is equal among species 
    δz_vec = range(0.05,stop=0.15,length=20) #Distance between species in resource space 
    c_vec = [0.1 0.15 0.2 0.3 0.4 0.5] #Intraspecific transmission rate coefficint
    
    #maximum intraspecific basic reproduction number
    R₀_aa_max = 2.0 
    R₀_bb_max = 1.5 

    γ = 0.1 #Recovery rate
    z_start = 0.18
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2
    Nₚ = 1500    

    #Parameters (Maximum intraspecific transmission rate)
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
        
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
        
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    t_end = 8000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]      
    Threads.@threads for i in eachindex(c_vec)
        for j in eachindex(δz_vec)
            δz =  δz_vec[j]
            c = c_vec[i]
            file_name = "output_c_$(i)_dz_$(j)"
            
            β_ab_max = (β_aa_max+β_bb_max)/2*c
            μ_b = μ_a+δz        
            z_end = 0.2+δz+0.05

            option = z_start,z_end,Nₚ,u₀,tspan

            syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude)

            strats = singular_strategies(syspar,option)
            
            JLD2.@save "./output/AD_2_Species/parameter_influence/model_check/"*file_name*".jld2" strats syspar
        end
    end

    c_len = length(c_vec)
    dz_len = length(δz_vec)

    JLD2.@save "./output/AD_2_Species/parameter_influence/model_check/c_dz_info.jld2" c_len dz_len
end

compute() 