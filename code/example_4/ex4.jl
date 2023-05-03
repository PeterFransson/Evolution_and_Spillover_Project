#Test: stocastic simulation using Gillespie algorithm [1] 
#
#References
#[1] Gillespie, D. T. (1977). Exact stochastic simulation of
#   coupled chemical reactions. The journal of physical chemistry,
#   81(25), 2340-2361.

function SI_Gillespie(S,I,t,p)
    c,γ,K,N = p

    rᵣ = I*c
    rᵢ = S*γ*K*I/N

    θ = 1/(rᵢ+rᵣ)    

    d = Exponential(θ) 
    Δt = rand(d)

    state = sample(['i','r'],Weights([rᵢ*θ,rᵣ*θ]))

    if state=='i'
        return S-1,I+1,t+Δt 
    elseif state=='r'
        return S+1,I-1,t+Δt 
    else
        error("State error: sate either i or r")
    end
end

function run_Gillespie(S₀,I₀,tspan,p)    
    t_end = tspan[2]
    S = S₀
    I = I₀
    t = tspan[1]
    t_vec = [t]
    S_vec = [S]
    I_vec = [I₀]

    while t<t_end
        S,I,t = SI_Gillespie(S,I,t,p)        

        push!(S_vec,S)
        push!(I_vec,I)
        push!(t_vec,t)        
    end
    return (S_vec,I_vec,t_vec)
end

function run_ex4()
    c = 0.1 #Recovery rate
    γ = 0.6 #Probability of infection    
    K = 0.45 #Contact rate   
    N = 10^5 #Population size 
    
    I₀ = 10
    S₀ = N-I₀   
    t₀ = 0.0
    
    p = (c,γ,K,N)
    @time S_vec,I_vec,t_vec = run_Gillespie(S₀,I₀,(t₀,80),p) 

    plt = plot(t_vec,S_vec,ylabel="time",legend=false)
    plot!(plt,t_vec,I_vec)

    for i = 1:3
        @time S_vec,I_vec,t_vec = run_Gillespie(S₀,I₀,(t₀,80),p)    
        plot!(plt,t_vec,S_vec)  
        plot!(plt,t_vec,I_vec)
    end 
    plot(plt)
end

run_ex4()