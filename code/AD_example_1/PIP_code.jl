function draw_PIP()
    Nₚ = 200
    z_vec = range(0.0,stop=1.0,length=Nₚ) 

    r_m_mat = zeros(Nₚ,Nₚ)

    #Parameters
    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3 

    c = 0.1    

    #Initial values    
    N_a = 10^3  
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 0
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
    t_end = 900
    tspan = (t_start,t_end)

    u₀ = [S₀_a S₀_b;I_a₀ I_b₀]

    mat_row_idx = collect(Nₚ:-1:1) 
    for col_idx in 1:1:Nₚ #Resident strain
        resident_strain = z_vec[col_idx]
        p = c,γ_a,γ_b,resident_strain,K_aa,K_bb,K_ab,N_a,N_b
        eq_point = find_eq(u₀,tspan,p)

        S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
        I_a_eq,I_b_eq = eq_point.I_a,eq_point.I_b

        I = I_a_eq+I_b_eq
        q_a = I_a_eq/I
        q_b = I_b_eq/I

        K_N_mean_a =  K_aa/N_a*q_a+K_ab/N_b*q_b
        K_N_mean_b = K_ab/N_a*q_a+K_bb/N_b*q_b 
    
        r_m(z) = γ_a(z)*S_a_eq*K_N_mean_a+γ_b(z)*S_b_eq*K_N_mean_b-c 

        for i in 1:1:Nₚ #Mutant strain strain
            row_idx = mat_row_idx[i]
            mutant_strain = z_vec[i]            
            r_m(mutant_strain)<0.0||(r_m_mat[i,col_idx]=1)         
        end        
    end

    heatmap(z_vec,z_vec,r_m_mat)
    plot!([0.0,1.0],[0.0,1.0])
    plot!([μ_a,μ_a],[0.0,1.0])
    plot!([μ_b,μ_b],[0.0,1.0])
end

draw_PIP()