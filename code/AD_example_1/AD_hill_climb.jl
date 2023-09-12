function hill_climb()
    N_steps = 1000 #Number of trait substitutions
    Nₚ = 1000
    z_start = 0.1
    z_end = 0.4
    z_vec = range(z_start,stop=z_end,length=Nₚ) 

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
    μ_b,σ²_b,amplitude_b = 0.32,0.0025,0.6

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

    steps = collect(1:N_steps)
    resident_traits = zeros(size(steps))    
    resident_traits[1] = strategy
    for i in 1:N_steps-1
        p = c,γ_a,γ_b,resident_traits[i],K_aa,K_bb,K_ab,N_a,N_b
        @show R0_val = calc_R0(p)
        
        (strategy_new,u₀_new) =find_new_strategy(u₀,tspan,p,γ_prim)
        
        resident_traits[i+1] = strategy_new 
        u₀ = u₀_new        
    end

    plot(steps,resident_traits)
end
hill_climb()