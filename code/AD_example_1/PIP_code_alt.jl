function calc_R0(p)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)
    K_L = [γ_aa*K_aa/c γ_aa*K_ab*N_a/(N_b*c);γ_bb*K_ab*N_b/(N_a*c) γ_bb*K_bb/c]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end

function calc_R0_m(p,S_a_eq,S_b_eq,strategy_m)
    c,γ_a,γ_b,strategy_r,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_aa = γ_a(strategy_m)
    γ_bb = γ_b(strategy_m)
    K_L = [γ_aa*K_aa*S_a_eq/(c*N_a) γ_aa*K_ab*S_a_eq/(N_b*c);γ_bb*K_ab*S_b_eq/(N_a*c) γ_bb*K_bb*S_b_eq/(c*N_b)]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end

function draw_PIP()
    Nₚ = 1000
    z_start = 0.1
    z_end = 0.4
    z_vec = range(z_start,stop=z_end,length=Nₚ) 

    r_m_mat = zeros(Nₚ,Nₚ)

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
    μ_b,σ²_b,amplitude_b = 0.3,0.0025,0.6

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b)
    γ_a_prim(z) = γ_fun_prim(z,μ_a,σ²_a,amplitude_a)
    γ_b_prim(z) = γ_fun_prim(z,μ_b,σ²_b,amplitude_b)

    γ_prim = γ_a_prim,γ_b_prim

    strategy = 0.2

    p = c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b

    t_start = 0
    t_end = 1500
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
                if calc_R0_m(p,S_a_eq,S_b_eq,mutant_strain)>1.0
                    r_m_mat[i,col_idx]=1 
                end                  
            else
                r_m_mat[i,col_idx]=0.5
            end                 
        end        
    end

    heatmap(z_vec,z_vec,r_m_mat)
    plot!([z_start,z_end],[z_start,z_end])
    plot!([μ_a,μ_a],[z_start,z_end])
    plot!([μ_b,μ_b],[z_start,z_end],legend=false)
    savefig("./fig/AD_example_1/pip_alt.svg")    
end

draw_PIP()