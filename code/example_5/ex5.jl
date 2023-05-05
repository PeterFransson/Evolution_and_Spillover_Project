#Test: stocastic one species evolutionary SI-model simulation using Gillespie algorithm [1] 
#

mutable struct Infected
    number::Integer
    strategy::Real
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
τ_fun(r,rand_num) = -log(rand_num)/r

function SI_Gillespie!(S,I::Vector{Infected},t,p)
    c,μₘ,σₘ,γ,K,N,t_end = p
    if isempty(I)==false
        n_I = length(I)

        r_num = rand(n_I,3)
        Δt = Inf
        idx = 1
        state = 1

        for i in 1:n_I
            β = K*γ(I[i].strategy)
            I_z = I[i].number
            rᵣ = I_z*c
            rₘ = I_z*μₘ
            rᵢ = S*β*I_z/N
            τ = [τ_fun(rᵣ,r_num[i,1]),τ_fun(rₘ,r_num[i,2]),τ_fun(rᵢ,r_num[i,3])]
            state_temp = argmin(τ)
            if τ[state_temp]<Δt
                Δt = τ[state_temp]
                idx = i
                state = state_temp
            end
        end        

        if state==1 #recover
            I_z = I[idx].number
            if I_z-1>0
                I[idx].number = I_z-1
            else
                popat!(I,idx)
            end
            return S+1,t+Δt 
        elseif state==2 #mutation     
            I_z = I[idx].number       
            resident_strategy = I[idx].strategy
            dₘ = Normal(resident_strategy,σₘ)
            mutant_strategy = min(max(rand(dₘ),0.0),1.0)
            push!(I,Infected(1,mutant_strategy))
            if I_z-1>0
                I[idx].number = I_z-1
            else
                popat!(I,idx)
            end
            return S,t+Δt 
        elseif state==3 #infection
            I[idx].number += 1
            return S-1,t+Δt 
        else
            error("State error: sate either i or r")
        end
    else
        return S,t_end
    end
end

function run_Gillespie!(S₀,I,t₀,p,n_samples)    
    t_end = p[7]
    S = S₀    
    t = t₀
    t_vec = [t]
    S_vec = [S]
    I_vec = [sum([inf.number for inf in I])]
    n_variants = [length(I)]
    mean_strategy = [mean([inf.strategy for inf in I])]

    Δt_sample = (t_end-t₀)/n_samples
    sample_nr = 1

    while t<t_end
        S,t = SI_Gillespie!(S,I,t,p)    

        if t>t₀+sample_nr*Δt_sample
            push!(S_vec,S)            
            push!(t_vec,t) 
            push!(n_variants,length(I))
            if  length(I)>0
                push!(mean_strategy,mean([infected.strategy for infected in I]))
                push!(I_vec,sum([infected.number for infected in I]))
            else
                push!(mean_strategy,0.0)
                push!(I_vec,0)
            end
            sample_nr += 1
        end       
    end
    return (S_vec,I_vec,t_vec,n_variants,mean_strategy)
end

function run_ex5()
    Nₚ = 200
    z_vec = range(0.0,stop=1.0,length=Nₚ+1)

    μ,σ²,amplitude = 0.4,0.0025,0.6 #0.3

    γ(z) = γ_fun(z,μ,σ²,amplitude)
    
    
    #Plot strategy
    plot(z_vec,γ.(z_vec),label="Species A",
    xlabel="Strategy",ylabel="γ")    

    #Parameters
    K = 0.425
    c = 0.1
    μₘ = 0.05
    
    z_start = 0.33    

    #Initial states
    N = 2*10^3     
    I = [Infected(2,z_start)]  
    S₀ = N-2  
    σₘ = 0.0158

    I_eq =  (1-c/(γ(μ)*K))*N
    
    @show γ.(z_start)*K/c

    #Plot R₀
    plot(z_vec,γ.(z_vec)*K/c,xlabel="Strategy",ylabel="R₀") 

    t₀ = 0.0
    t_end = 200.0

    p = c,μₘ,σₘ,γ,K,N,t_end

    n_samples = 500

    S_vec,I_vec,t_vec,n_variants,mean_strategy = run_Gillespie!(S₀,I,t₀,p,n_samples) 

    pl1 = plot(t_vec,mean_strategy)
    pl1 = plot!([t₀,t_end],[μ,μ])
    pl2 = plot(t_vec,n_variants)
    pl2 = plot!([t₀,t_end],[I_eq,I_eq])

    pl3 = plot(t_vec,S_vec)
    pl3 = plot!(t_vec,I_vec)
    pl3 = plot!([t₀,t_end],[I_eq,I_eq])

    l = @layout [a b;c]

    plot(pl1,pl2,pl3,layout=l,legends=false)
end

run_ex5()