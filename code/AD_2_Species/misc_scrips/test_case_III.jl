function run_case_III_test()
    Nₚ = 2000
    z_start = 0.18#0.1
    z_end = 0.4#0.5

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.35#0.35
    
    N_tot = 2*10^3 
    #r_N = 0.1
    r_N = 0.5
    N_a = N_tot*r_N#10^3  
    N_b = N_tot*(1-r_N)#10^3

    syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀  

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    option = z_start,z_end,Nₚ,u₀,tspan    

    strats = singular_strategies(syspar,option)

    z = strats[1].strategy

    @show R₀_a = τ_fun(μ_a,z,σ²,amplitude)*c_aa/γ
    @show R₀_b = τ_fun(μ_b,z,σ²,amplitude)*c_bb/γ

    p = (z,syspar)
    find_eq(u₀,tspan,p)
end

run_case_III_test()