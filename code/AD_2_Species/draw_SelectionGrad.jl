function calc_selectgrad(strain,p_in,option)
    u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b = p_in

    p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
    eq_point = find_eq(u₀,tspan,p_eq)

    p_grad = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain)
    S = (eq_point.S_a,eq_point.S_b)
    dλ_max = calculate_select_grad(p_grad,S)

    return dλ_max
end
function calc_d_selectgrad(strain,p_in,option)
    u₀,tspan = option    
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b = p_in

    p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
    eq_point = find_eq(u₀,tspan,p_eq)

    p_curv = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain),τ_d_prime_a(strain),τ_d_prime_b(strain)
    S = (eq_point.S_a,eq_point.S_b)
    ddλ_max = calculate_sec_inv_fit(p_curv,S)

    return ddλ_max
end

function calc_selectgrad(p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b = p_in

    strategies = range(z_start,stop=z_end,length=Nₚ) 
    select_grad = zeros(Nₚ)

    for i = 1:Nₚ
        strain = strategies[i]
        p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
        eq_point = find_eq(u₀,tspan,p_eq)

        p_grad = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain)
        S = (eq_point.S_a,eq_point.S_b)
        dλ_max = calculate_select_grad(p_grad,S)

        select_grad[i] = dλ_max
    end

    return (select_grad,strategies)
end

function plot_def()
    Nₚ = 1000
    z_start = 0.18
    z_end = 0.36

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
    τ_d_prime_a(z) = τ_d_prime_fun(z,μ_a,σ²,amplitude)
    τ_d_prime_b(z) = τ_d_prime_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]     

    option = z_start,z_end,Nₚ,u₀,tspan
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b

    select_grad,strategies = calc_selectgrad(p_in,option)

    option_bisec = u₀,tspan
    z_0 = bisection(x->calc_selectgrad(x,p_in,option_bisec), z_start,z_end)
    @show z_0
    @show calc_selectgrad(z_0,p_in,option_bisec)
    p_d_selectgrad = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b
    @show ddλ_max_0 = calc_d_selectgrad(z_0,p_d_selectgrad,option_bisec)

    plot([z_start,z_end],[0,0])#,ylim=[-0.2,0.2])
    plot!(strategies,select_grad,ylabel="Sᵣ(m=r)′",xlabel="Strain",legends=false)
end

plot_def()