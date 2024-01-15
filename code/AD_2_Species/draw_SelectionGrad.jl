function draw_gif()
    δz_vec = range(0.05,stop=0.15,length=100)

    γ = 0.1 
    z_start = 0.18
    σ²,amplitude = 0.0025,0.6
    μ_a = 0.2
    Nₚ = 1500    

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

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    
    
    select_grad_gif = @animate for δz in δz_vec
        
        μ_b = 0.2+δz        
        z_end = 0.2+δz+0.05

        option = z_start,z_end,Nₚ,u₀,tspan
        syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)
        
        strats = singular_strategies(syspar,option)
        
        strategies = range(z_start,stop=z_end,length=200)
        select_grad = find_selectgrad.(strategies,Ref(syspar),Ref((u₀,tspan)))

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

draw_gif()