function test()
    Nₚ = 1000
    z_start = 0.1
    z_end = 0.4

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.35    
    
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)
    τ_prime_a(z) = τ_prime_fun(z,μ_a,σ²,amplitude)
    τ_prime_b(z) = τ_prime_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]  
    
    r = 0.3475

    option = z_start,z_end,Nₚ,u₀,tspan     
    p_sys = τ_a(r),τ_b(r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    eq = find_eq(u₀,tspan,p_sys)


    p_grad = τ_a(r),τ_b(r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a(r),τ_prime_b(r)
    S = (eq.S_a,eq.S_b)
    calculate_select_grad(p_grad,S)
end

test()