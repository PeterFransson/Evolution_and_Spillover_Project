#Test: stocastic SI-model simulation using Gillespie algorithm [1] 
#
#References
#[1] Gillespie, D. T. (1977). Exact stochastic simulation of
#   coupled chemical reactions. The journal of physical chemistry,
#   81(25), 2340-2361.

function SI!(du,u,p,t)
    c,β,N,t_end = p
    S,I = u
    du[1] = -S*β*I/N+I*c
    du[2] = -du[1]
end 

function SI_Gillespie(S,I,t,p)
    c,β,N,t_end = p
    if I>0
        rᵣ = I*c
        rᵢ = S*β*I/N

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
    else
        return S,I,t_end
    end
end

function run_Gillespie(S₀,I₀,t₀,p,n_samples)    
    t_end = p[4]
    S = S₀
    I = I₀
    t = t₀
    t_vec = [t]
    S_vec = [S]
    I_vec = [I₀]

    Δt_sample = (t_end-t₀)/n_samples
    sample_nr = 1

    while t<t_end
        S,I,t = SI_Gillespie(S,I,t,p)    

        if t>t₀+sample_nr*Δt_sample
            push!(S_vec,S)
            push!(I_vec,I)
            push!(t_vec,t) 
            sample_nr += 1
        end       
    end
    return (S_vec,I_vec,t_vec)
end

function run_ex4()
    c = 0.1 #Recovery rate
    β = 0.09 #Transmission rate    
    N = 10^3 #Population size 
    
    I₀ = 10
    S₀ = N-I₀   
    t₀ = 0.0
    t_end = 300.0
    
    p = (c,β,N,t_end)
    n_sample = 200
    @time S_vec,I_vec,t_vec = run_Gillespie(S₀,I₀,t₀,p,n_sample)     
            
    plt = plot(t_vec,S_vec,ylabel="time",legend=false)
    plot!(plt,t_vec,I_vec)  

    for i in 1:100
        @time S_vec,I_vec,t_vec = run_Gillespie(S₀,I₀,t₀,p,n_sample)     
        
        plot!(plt,t_vec,S_vec)
        plot!(plt,t_vec,I_vec)  
    end

    u₀ = [S₀,I₀]    
    prob = ODEProblem(SI!,u₀,(t₀,t_end),p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE
    t_vec = range(t₀,stop=t_end,length=100)
    S_vec = first.(sol.(t_vec))
    I_vec = last.(sol.(t_vec))
    plot!(plt,t_vec,S_vec,seriestype=:scatter)  
    plot!(plt,t_vec,I_vec,seriestype=:scatter)  

    plot(plt)    
end


run_ex4()