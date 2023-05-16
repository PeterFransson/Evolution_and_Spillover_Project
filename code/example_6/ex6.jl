#Test: stocastic two species evolutionary SI-model simulation using Gillespie algorithm [1] 
#

mutable struct Infected
    number::Vector{Integer}  
    strategy::Real
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
τ_fun(r,rand_num) = -log(rand_num)/r

function SI_Gillespie!(S::Vector{U},
    I::Vector{Infected},
    n_species::T,
    t,
    p) where {U<:Integer,T<:Integer}

    c,μₘ,σₘ,γ,K,N,t_end = p
    if isempty(I)==false
        n_pathogen_strains = length(I)        

        r_num = rand(n_pathogen_strains,n_species,3)
        Δt = Inf
        pathogen_strains = 1
        species = 1
        state = 1   

        for i in 1:n_pathogen_strains #Pathogen varriant loop
            I_z = I[i].number
            strategy_z = I[i].strategy
            for j in 1:n_species #species loop                       
                 
                rᵣ = I_z[j]*c
                rₘ = I_z[j]*μₘ
                rᵢ = S[j]*γ[j](strategy_z)*sum(K[j,:].*I_z./N)

                τ = [τ_fun(rᵣ,r_num[i,j,1]),τ_fun(rₘ,r_num[i,j,2]),τ_fun(rᵢ,r_num[i,j,3])]
                state_temp = argmin(τ)
                if τ[state_temp]<Δt
                    Δt = τ[state_temp]
                    pathogen_strains = i
                    state = state_temp
                    species = j
                end                
            end       
        end

        if state==1 #recover
            I_z = I[pathogen_strains].number
            if sum(I_z)-1>0
                I[pathogen_strains].number[species] -= 1
            else
                popat!(I,pathogen_strains)
            end
            S[species]+=1
            return t+Δt 
        elseif state==2 #mutation     
            I_z = I[pathogen_strains].number      
            resident_strategy = I[pathogen_strains].strategy
            dₘ = Normal(resident_strategy,σₘ)
            mutant_strategy = min(max(rand(dₘ),0.0),1.0)
            I_z_mutant = zeros(eltype(I_z),n_species)
            I_z_mutant[species]+=1
            push!(I,Infected(I_z_mutant,mutant_strategy))
            if sum(I_z)-1>0
                I[pathogen_strains].number[species] -= 1
            else
                popat!(I,pathogen_strains)
            end
            return t+Δt 
        elseif state==3 #infection
            I[pathogen_strains].number[species] += 1
            S[species]-=1
            return t+Δt 
        else
            error("State error: sate either i, m, or r")
        end
    else
        return t_end
    end
end

function run_Gillespie!(S::Vector{T},
    I::Vector{Infected},
    t₀,
    p,
    n_samples::Integer) where {T<:Integer}

    t_end = last(p)        
    t = t₀
    t_vec = [t]
    S_vec = typeof(S)[]
    push!(S_vec,S) ##

    I_vec = typeof(I)[] 
    push!(I_vec,copy(I))     
    n_variants = [length(I)]
    mean_strategy = [mean([strain.strategy for strain in I])]
    
    Δt_sample = (t_end-t₀)/n_samples
    sample_nr = 1
    n_species = length(S)
    while t<t_end
        t = SI_Gillespie!(S,I,n_species,t,p)    

        if t>t₀+sample_nr*Δt_sample
            println("$(sample_nr)/$(n_samples)")
            push!(S_vec,S) ##          
            push!(t_vec,t) 
            push!(n_variants,length(I))
            if  length(I)>0
                push!(mean_strategy,mean([infected.strategy for infected in I]))
                #push!(I_vec,sum([infected.number for infected in I]))
                push!(I_vec,copy(I))                  
            else
                push!(mean_strategy,0.0)
                push!(I_vec,Infected[])
            end
            sample_nr += 1            
        end       
    end    
    return (S_vec,I_vec,t_vec,n_variants,mean_strategy)
end

function draw_strain_distribution(t_vec,I_vec,img_name)
        
    strain_distribution = @animate for (time,strain_dist) in zip(t_vec,I_vec)      
        
        I_strain = [sum(strain.number) for strain in strain_dist]
        strategy_strain = [strain.strategy for strain in strain_dist]
        I_tot = sum(I_strain)
        q_strain = I_strain/I_tot
        
        plot(strategy_strain,
        q_strain,
        seriestype=:scatter,
        legends=false,
        ylims=(0,0.1),xlim=(0,1))        
    end
    gif(strain_distribution, "./fig/ex6/"*img_name*".gif", fps = 15)      
end

function draw_strain_evolution_plot(t_vec,I_vec,t₀,t_end,img_name)
    plt = plot(ylims=(0.0,1.0),xlims=(t₀,t_end),legends=false)
    for (time,strain_dist) in zip(t_vec,I_vec) 
        strategies = [strain.strategy for strain in strain_dist]
        time_vec = [time for strain in strain_dist]
        plot!(plt,time_vec,strategies,seriestype=:scatter)              
    end    
    savefig(plt,"./fig/ex6/"*img_name*".svg")
end

function run_ex6()
    Nₚ = 200
    z_vec = range(0.0,stop=1.0,length=Nₚ+1)

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.6
    μ_b,σ²_b,amplitude_b = 0.4,0.0025,0.6 #0.3

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b) 
   
    #Plot strategy
    plot(z_vec,γ_a.(z_vec),label="Species A",
    xlabel="Strategy",ylabel="γ") 
    plot(z_vec,γ_b.(z_vec),label="Species B")   

    #Parameters
    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3 
    c = 0.1
    μₘ = 0.05
    
    z_start = 0.2    

    #Initial states
    N_a = 10^3     
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 
    σₘ = 0.05#0.0158

    S = [S₀_a,S₀_b] 
    I = [Infected([I_a₀,I_b₀],z_start)]  

    @show I_a_eq =  (1-c/(γ_a(μ_a)*K_aa))*N_a
    
    @show γ_a.(z_start)*K_aa/c
    @show γ_b.(z_start)*K_bb/c
    
    #Plot R₀
    plot(z_vec,γ_a.(z_vec)*K_aa/c,xlabel="Strategy",ylabel="R₀") 
    plot!(z_vec,γ_b.(z_vec)*K_bb/c) 

    
    t₀ = 0.0
    t_end = 1500.0

    K = [K_aa K_ab;K_ab K_bb]
    γ = (γ_a,γ_b)
    N = [N_a,N_b]

    p = c,μₘ,σₘ,γ,K,N,t_end

    n_samples = 50

    S_vec,I_vec,t_vec,n_variants,mean_strategy = run_Gillespie!(S,I,t₀,p,n_samples) 
    
    draw_strain_distribution(t_vec,I_vec,"strain_evo")
    draw_strain_evolution_plot(t_vec,I_vec,t₀,t_end,"trait_vs_time")

    pl1 = plot(t_vec,mean_strategy)
    pl1 = plot!([t₀,t_end],[μ_a,μ_a])
    pl1 = plot!([t₀,t_end],[μ_b,μ_b])
    pl2 = plot(t_vec,n_variants)    
    
    plot(pl1,pl2,layout=(2,1),legneds=false)
    
    #=
    pl3 = plot(t_vec,S_vec)
    pl3 = plot!(t_vec,I_vec)
    pl3 = plot!([t₀,t_end],[I_eq,I_eq])

    l = @layout [a b;c]

    plot(pl1,pl2,pl3,layout=l,legends=false)
    =#
end

run_ex6()