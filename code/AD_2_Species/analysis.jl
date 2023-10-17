function test()
    Nₚ = 1000
    z_start = 0.1
    z_end = 0.5

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.32    
    
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 6000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    option = z_start,z_end,Nₚ,u₀,tspan 

    #p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b
    strategy = 0.32
    p_strategy = τ_a(strategy),τ_b(strategy),c_aa,c_bb,c_ab,γ,γ,N_a,N_b

    R = system_matrix(p_strategy,(S₀_a,S₀_b))
    @show eigvals(R)
    
    eq_point = find_eq(u₀,tspan,p_strategy)

    S_eq = (eq_point.S_a,eq_point.S_b)
    R = system_matrix(p_strategy,S_eq)
    @show eigvals(R)

    strategy_new = strategy-0.00001
    p_strategy = τ_a(strategy_new),τ_b(strategy_new),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    R = system_matrix(p_strategy,S_eq)
    @show eigvals(R)
end

test()