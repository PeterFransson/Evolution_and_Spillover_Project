function plot_def()
    folder_name = "./fig/AD_2_Species/Conv/between_case_b_c_2"
    name = "between_case_b_c_2"
    isdir(folder_name)|| mkdir(folder_name)

    Nₚ = 1500
    z_start = 0.18
    z_end = 0.33
    Δz = (z_end-z_start)/(Nₚ-1)

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.28
    
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
    p_strats = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b

    select_grad,strategies = calc_selectgrad(p_in,option)

    strats = singular_strategies(p_strats,option)
    @show strats 

    #=
    option_bisec = u₀,tspan
    
    z0s = find_singular_strategies(p_in,option)
    @show z0s
    conv_vecs = check_convergence.(z0s,Ref(p_in),Ref(option_bisec);h=Δz)
    @show conv_vecs
    ddλ_max_0s = calc_d_selectgrad.(z0s,Ref(p_d_selectgrad),Ref(option_bisec))
    @show ddλ_max_0s
    =#

    plot([z_start,z_end],[0,0])#,ylim=[-1.0,1.0])
    plot!(strategies,select_grad,ylabel="Sᵣ(m=r)′",xlabel="Strain",legends=false)
    #savefig(folder_name*"/"*name*".svg")
    #savefig("./fig/PIP/2_species/between_case_b_c/singular_strategies.svg")
end

function get_selct_grad_def(δz::Real)
    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a = 0.2
    μ_b = 0.2+δz

    Nₚ = 1500
    z_start = 0.18
    z_end = μ_b+0.05
    Δz = (z_end-z_start)/(Nₚ-1)

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3
    
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
    p_strats = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b

    select_grad,strategies = calc_selectgrad(p_in,option)
    
    strats = singular_strategies(p_strats,option)
    
    return (select_grad,strategies,strats)
end

function draw_gif()
    δz_vec = range(0.05,stop=0.15,length=100)
    
    select_grad_gif = @animate for δz in δz_vec
        μ_a = 0.2
        μ_b = 0.2+δz
        z_start = 0.18
        z_end = 0.2+δz+0.05
        select_grad,strategies,strats = get_selct_grad_def(δz)
        min_grad = minimum(select_grad)
        max_grad = maximum(select_grad)

        plot([z_start,z_end],[0,0])
        plot!([μ_a,μ_a],0.05*[min_grad,max_grad],color=:blue)
        plot!([μ_b,μ_b],0.05*[min_grad,max_grad],color=:red)
        plot!(strategies,select_grad,ylabel="Sᵣ(m=r)′",xlabel="Strain",legends=false)
        for strat in strats 
            if strat.evo_stable
                scatter!([strat.strategy],[0],mc=:green, ms=5)
            else
                scatter!([strat.strategy],[0],mc=:red, ms=5)
            end
        end
    end
    gif(select_grad_gif, "./fig/AD_2_Species/select_grad.gif", fps = 5) 
end

#plot_def()
draw_gif()