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
    τ_prim_a(z) = τ_prim_fun(z,μ_a,σ²,amplitude)
    τ_prim_b(z) = τ_prim_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    option = z_start,z_end,Nₚ,u₀,tspan 
    
    #=
    z_r = 0.25
    p_system = τ_a(z_r),τ_b(z_r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    eq = find_eq(u₀,tspan,p_system)
    S = (eq.S_a,eq.S_b)

    R = system_matrix(p_system,S)
    @show R_eval = maximum(eigvals(R))

    Δz = 0.00001
    z_m = z_r+Δz 
    p_system = τ_a(z_m),τ_b(z_m),c_aa,c_bb,c_ab,γ,γ,N_a,N_b
    R = system_matrix(p_system,S)
    @show R_eval_new = maximum(eigvals(R))

    @show dλ_max_num = (R_eval_new-R_eval)/Δz

    p_grad = τ_a(z_r),τ_b(z_r),c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prim_a(z_r),τ_prim_b(z_r)
    calculate_select_grad(p_grad,S) 
    =# 

    option = z_start,z_end,Nₚ,u₀,tspan 
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b 

    coex_region = create_coex_region(p_in,option)   
    z_vec = range(z_start,stop=z_end,length=Nₚ) 
    heatmap(z_vec,z_vec,coex_region)

    u₀_2_strain = [N_a-2 1 1;N_b-2 1 1]   

    option = z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b,τ_prim_a,τ_prim_b 
    tep = create_TEP(p_in,option,coex_region;ϵ=1e-2)

    heatmap(z_vec,z_vec,tep) 
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    plot!([z_start,z_end],[μ_a,μ_a],c=:red)
    plot!([z_start,z_end],[μ_b,μ_b],c=:blue)
    plot!([μ_b,μ_b],[z_start,z_end],legend=false,c=:blue)
    savefig("./fig/AD_2_Species/tep.svg")    
end

test()