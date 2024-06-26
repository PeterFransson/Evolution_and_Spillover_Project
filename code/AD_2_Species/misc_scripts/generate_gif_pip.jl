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

function generate_gif_pip()
    R₀_ratio = 1.1#[1.0,1.5,2.0]
    c_ratio = 0.8#[0.8,0.6,0.2]
    Δᵣ_vec = collect(range(1.0,stop=2.5,length=50))   

    pip_gif = @animate for Δᵣ in Δᵣ_vec

        syspar = create_syspar(R₀_ratio,Δᵣ,c_ratio) 
        
        z_start = syspar.z_a-0.02
        z_start>0||error("z_start<0")
        z_end = syspar.z_b+0.02
        
        #Initial values
        Nₚ = 1000
        N_a = syspar.N_a    
        N_a>2||error("N_a<3")
        I_a₀ = 1  
        S₀_a = N_a-I_a₀ 
         
        N_b = syspar.N_b   
        N_b>2||error("N_b<3")
        I_b₀ = 1
        S₀_b = N_b-I_b₀  
    
        t_start = 0
        t_end = 7000
        tspan = (t_start,t_end)
    
        u₀ = [S₀_a I_a₀;S₀_b I_b₀]  
    
        option = z_start,z_end,Nₚ,u₀,tspan 

        strategies = range(z_start,stop=z_end,length=Nₚ)
        
        r_m_mat = create_PIP(syspar,option)

        z_a,z_b = syspar.z_a,syspar.z_b

        heatmap(strategies,strategies,r_m_mat)
        plot!([z_start,z_end],[z_start,z_end],c=:green)
        plot!([z_a,z_a],[z_start,z_end],c=:red)
        plot!([z_b,z_b],[z_start,z_end],xlims=[z_start,z_end],legend=false,c=:blue)
    end

    gif(pip_gif,"./fig/AD_2_Species/R0_slight_dif.gif", fps = 5) 
end

generate_gif_pip()