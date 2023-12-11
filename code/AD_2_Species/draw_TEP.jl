function draw_TEP()
    Nₚ = 1500
    z_start = 0.1
    z_end = 0.35

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.323    
    
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
    τ_d_prime_a(z) = τ_d_prime_fun(z,μ_a,σ²,amplitude)
    τ_d_prime_b(z) = τ_d_prime_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    
    u₀_2_strain = [N_a-2 1 1;N_b-2 1 1]  

    #option = z_start,z_end,Nₚ,u₀,tspan 
    
    #=
    z_r = 0.3
    p_system = τ_a(z_r),τ_b(z_r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    eq = find_eq(u₀,tspan,p_system)
    S = (eq.S_a,eq.S_b)

    R = system_matrix(p_system,S)
    @show R_eval = maximum(eigvals(R))

    Δz = 0.0000001
    z_m = z_r+Δz 
    p_system = τ_a(z_m),τ_b(z_m),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    R = system_matrix(p_system,S)
    @show R_eval_R = maximum(eigvals(R))
    
    z_m = z_r-Δz 
    p_system = τ_a(z_m),τ_b(z_m),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    R = system_matrix(p_system,S)
    @show R_eval_L = maximum(eigvals(R))

    @show dλ_max_num = (R_eval_R-R_eval)/Δz
    @show ddλ_max_num = (R_eval_R-2*R_eval+R_eval_L)/Δz^2

    p_grad = τ_a(z_r),τ_b(z_r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a(z_r),τ_prime_b(z_r)
    @show d_R = calculate_select_grad(p_grad,S)   

    z_m = z_r+Δz 
    p_grad = τ_a(z_m),τ_b(z_m),c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a(z_m),τ_prime_b(z_m)
    @show d_R_R = calculate_select_grad(p_grad,S)
    
    @show (d_R_R-d_R)/Δz
    
    p_sec = τ_a(z_r),τ_b(z_r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a(z_r),τ_prime_b(z_r),τ_d_prime_a(z_r),τ_d_prime_b(z_r)
    @show calculate_sec_inv_fit(p_sec,S) 
    =#
    
    option = z_start,z_end,Nₚ,u₀,tspan 
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b 

    coex_region = create_coex_region(p_in,option)   
    z_vec = range(z_start,stop=z_end,length=Nₚ) 
    heatmap(z_vec,z_vec,coex_region)     

    
    option = z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b 
    tep = create_TEP(p_in,option,coex_region;ϵ=1e-2)
    
    option = z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b,τ_prime_a,τ_prime_b

    #Check cone of invasions  
    @show calc_invasion_cone(0.174,0.348,p_in,option)          
    @show calc_invasion_cone(0.178,0.286,p_in,option)   
    @show calc_invasion_cone(0.229,0.340,p_in,option)    
    @show calc_invasion_cone(0.228,0.289,p_in,option)    
    
    heatmap(z_vec,z_vec,tep) 
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    plot!([z_start,z_end],[μ_a,μ_a],c=:red)
    plot!([z_start,z_end],[μ_b,μ_b],c=:blue)
    plot!([μ_b,μ_b],[z_start,z_end],legend=false,c=:blue)
    savefig("./fig/TEP/2_species/between_case_b_c/TEP.svg")  
     
end

draw_TEP()