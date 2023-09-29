function draw_PIP()
    Nₚ = 1000
    z_start = 0.1
    z_end = 0.5
    z_vec = range(z_start,stop=z_end,length=Nₚ) 

    r_m_mat = zeros(Nₚ,Nₚ)
    r_m_mat_alt = zeros(Nₚ,Nₚ)
    r_m_mat_alt_alt = zeros(Nₚ,Nₚ)

    #Parameters
    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3#0.3 

    c = 0.1    

    #Initial values    
    N_a = 10^3  
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 1
    S₀_b = N_b-I_b₀ 

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.6
    μ_b,σ²_b,amplitude_b = 0.35,0.0025,0.6

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b)
    γ_a_prim(z) = γ_fun_prim(z,μ_a,σ²_a,amplitude_a)
    γ_b_prim(z) = γ_fun_prim(z,μ_b,σ²_b,amplitude_b)

    γ_prim = γ_a_prim,γ_b_prim

    strategy = 0.2

    p = c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b

    t_start = 0
    t_end = 5500
    tspan = (t_start,t_end)

    u₀ = [S₀_a S₀_b;I_a₀ I_b₀]       

    for col_idx in 1:1:Nₚ #Resident strain
        resident_strain = z_vec[col_idx]
        p = c,γ_a,γ_b,resident_strain,K_aa,K_bb,K_ab,N_a,N_b
        @show R0_val = calc_R0(p)
        eq_point = find_eq(u₀,tspan,p)

        S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
        I_a_eq,I_b_eq = eq_point.I_a,eq_point.I_b

        I = I_a_eq+I_b_eq
        q_a = I_a_eq/I
        q_b = I_b_eq/I        

        K_N_mean_a(q_a,q_b) = K_aa/N_a*q_a+K_ab/N_b*q_b
        K_N_mean_b(q_a,q_b) = K_ab/N_a*q_a+K_bb/N_b*q_b             
        
        r_m(z,q_a,q_b) = γ_a(z)*S_a_eq*K_N_mean_a(q_a,q_b)+γ_b(z)*S_b_eq*K_N_mean_b(q_a,q_b)-c

        for i in 1:1:Nₚ #Mutant strain strain
            mutant_strain = z_vec[i]  
            
            if R0_val>1.0
                if γ_a(mutant_strain)>=γ_a(resident_strain)&&r_m(mutant_strain,1.0,0.0)>0.0
                    r_m_mat[i,col_idx]=1 
                end  

                if γ_b(mutant_strain)>=γ_b(resident_strain)&&r_m(mutant_strain,0.0,1.0)>0.0
                    r_m_mat[i,col_idx]=1
                end                 
                if calc_R0_m(p,S_a_eq,S_b_eq,mutant_strain)>1.0
                    r_m_mat_alt[i,col_idx]=1
                end
                if check_system_matrix(p,mutant_strain,(S_a_eq,S_b_eq))
                    r_m_mat_alt_alt[i,col_idx] = 1
                end
            else
                r_m_mat[i,col_idx]=0.5
                r_m_mat_alt[i,col_idx]=0.5
                r_m_mat_alt_alt[i,col_idx] = 0.5
            end                 
        end        
    end

    heatmap(z_vec,z_vec,r_m_mat)
    plot!([z_start,z_end],[z_start,z_end])
    plot!([μ_a,μ_a],[z_start,z_end])
    plot!([μ_b,μ_b],[z_start,z_end],legend=false)
    savefig("./fig/AD_example_1/pip.svg")    

    heatmap(z_vec,z_vec,r_m_mat_alt)
    plot!([z_start,z_end],[z_start,z_end])
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    z_crit_a = μ_a+sqrt(2*σ²_a*log(amplitude_a*K_aa/c))
    @show γ_a(z_crit_a)*K_aa/c
    plot!([z_crit_a,z_crit_a],[z_start,z_end],linestyle=:dash,c=:red)
    z_crit_b = μ_b-sqrt(2*σ²_b*log(amplitude_b*K_bb/c))
    @show γ_b(z_crit_b)*K_bb/c
    plot!([z_crit_b,z_crit_b],[z_start,z_end],linestyle=:dash,c=:blue)
    plot!([μ_b,μ_b],[z_start,z_end],legend=false,c=:blue)
    savefig("./fig/AD_example_1/pip_alt.svg") 

    heatmap(z_vec,z_vec,r_m_mat_alt_alt)
    plot!([z_start,z_end],[z_start,z_end])
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    plot!([μ_b,μ_b],[z_start,z_end],legend=false,c=:blue)
    savefig("./fig/AD_example_1/pip_alt_alt.svg") 

    p = c,γ_a,γ_b,0.2,K_aa,K_bb,K_ab,N_a,N_b
    draw_vector_field(p,1.0,N_a-1,20,1.0,N_b,20)    
end

draw_PIP()