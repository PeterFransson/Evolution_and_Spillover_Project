#Test: stocastic three species evolutionary SI-model simulation using Gillespie algorithm
#with simple within-host competition [updated algorithm] discrete

mutable struct Infected
    number::Vector{Integer}  
    strategy::Real
end

#Simpson's 1/3 rule
simpson13(f,a,b) = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6
#Discretize a continuous positive function f
function cont2disc(f,disc_samples)
    n = length(disc_samples)
    areas = simpson13.(Ref(f),disc_samples[1:n-1],disc_samples[2:n])
    area_tot = sum(areas)
    return areas/area_tot
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
τ_fun(r,rand_num) = -log(rand_num)/r

function find_weighted_idx(rnd::R,c_prob::Vector{R}) where {R<:Real}
    #c_prob is the cumulative probablity vector for all indices
    idx = findfirst(rnd.<c_prob)
    return idx
end

#Calculate strategy value of pathogen_idx∈[1,2,…,nₚ]
z_fun(pathogen_idx::Integer,nₚ::Integer) = (pathogen_idx-0.5)/nₚ

#Create mutation probabilities. Note pathogen strategy z∈[0,1]
function create_mutation_matrix(nₚ::T,σₘ::R) where {T<:Integer,R<:Real}
    m = zeros(nₚ,nₚ)
    c_m = zeros(nₚ,nₚ)
    disc_samples = [(i-1)/nₚ for i in 1:nₚ+1]
    for i in 1:nₚ
        μₘ = z_fun(i,nₚ) 
        f(x) = normal_d(x,μₘ,σₘ) 
        m_prob = cont2disc(f,disc_samples)
        c_m_prob = cumsum(m_prob)

        m[i,:] = m_prob
        c_m[i,:] = c_m_prob
    end
    return (m,c_m)
end 

function SI_Gillespie!(S::Vector{U},
    I::Vector{Infected},
    n_species::T,
    t,
    p) where {U<:Integer,T<:Integer}

    c,μₘ,σₘ,γ,K,N,c_m,t_end = p
    
        n_pathogen_strains = length(I)        

        r_num = rand(n_pathogen_strains,n_species,3)
        Δt = Inf
        pathogen_strains = 0
        species = 0
        state = 0       

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
        I[pathogen_strains].number[species] -= 1            
        S[species]+=1
        return t+Δt 
    elseif state==2 #mutation               
        resident_strategy = I[pathogen_strains].strategy

        r_rnd = rand()            
        pathogen_strains_mutant = find_weighted_idx(r_rnd,c_m[pathogen_strains,:]) 
        mutant_strategy = I[pathogen_strains_mutant].strategy

        if γ[species](mutant_strategy)>γ[species](resident_strategy)
            I[pathogen_strains_mutant].number[species]+=1
            I[pathogen_strains].number[species] -= 1                
        end
        return t+Δt 
    elseif state==3 #infection
        I[pathogen_strains].number[species] += 1
        S[species]-=1            

        return t+Δt 
    else
        println("State error or extinction: sate either i, m, or r")
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
    push!(S_vec,copy(S))

    I_vec = typeof(I)[] 
    push!(I_vec,deepcopy(I))     
    
    Δt_sample = (t_end-t₀)/n_samples
    sample_nr = 1
    n_species = length(S)

    while t<t_end
        t = SI_Gillespie!(S,I,n_species,t,p)    

        if t>t₀+sample_nr*Δt_sample
            println("$(sample_nr)/$(n_samples)")
            push!(S_vec,copy(S))         
            push!(t_vec,t) 
            push!(I_vec,deepcopy(I))          
            sample_nr += 1            
        end       
    end 
        
    return (S_vec,I_vec,t_vec)
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
        xlim=(0,1))        
    end
    gif(strain_distribution, "./fig/ex9/"*img_name*".gif", fps = 15)      
end

function draw_strain_evolution_plot(t_vec,I_vec,t₀,t_end,img_name)
    plt = plot(ylims=(0.0,1.0),xlims=(t₀,t_end),legends=false)
    for (time,strain_dist) in zip(t_vec,I_vec)
        
        strategies = []
        for strain in strain_dist
            if sum(strain.number)>0
                push!(strategies,strain.strategy)
            end
        end

        if isempty(strategies) == false
            time_vec = [time for i in strategies]
            plot!(plt,time_vec,strategies,seriestype=:scatter) 
        end             
    end    
    savefig(plt,"./fig/ex9/"*img_name*".svg")
end

function draw_compartments(t_vec,I_vec,S_vec,img_name)    
    n_species = length(S_vec[1])
    n_time = length(t_vec)
    I_plot_mat = zeros(n_time,n_species)
    S_plot_mat = zeros(n_time,n_species)
    for t_idx in 1:n_time
        n_strains = length(I_vec[t_idx])
        for strain_idx in 1:n_strains            
            I_plot_mat[t_idx,:]+=I_vec[t_idx][strain_idx].number
        end
        S_plot_mat[t_idx,:]+=S_vec[t_idx]
    end    
    plt = plot(xlabel="Time",legend=false)
    for species_idx in 1:n_species
        plot!(plt,t_vec,S_plot_mat[:,species_idx])
        plot!(plt,t_vec,I_plot_mat[:,species_idx])
    end
    savefig(plt,"./fig/ex9/"*img_name*".svg")
end

function run_ex9()    
    Nₚ = 100 #Number of strains 
    z_vec = z_fun.(1:Nₚ,Ref(Nₚ))

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.6
    μ_b,σ²_b,amplitude_b = 0.3,0.0025,0.6 #0.3
    μ_c,σ²_c,amplitude_c = 0.5,0.0025,0.6

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b) 
    γ_c(z) = γ_fun(z,μ_c,σ²_c,amplitude_c)
   
    #Plot strategy
    plot(z_vec,γ_a.(z_vec),label="Species A",
    xlabel="Strategy",ylabel="γ") 
    plot!(z_vec,γ_b.(z_vec),label="Species B")   
    plot!(z_vec,γ_c.(z_vec),label="Species C")    
    
    #Parameters
    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3 
    K_cc = 0.425
    K_ac = 0.3 
    K_bc= 0.3 

    c = 0.1
    μₘ = 0.05

    @show z_start_idx = floor(Int64,0.2*Nₚ)
    @show z_start = z_fun(z_start_idx,Nₚ) 

    #Initial states
    N_a = 10^3     
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 
    N_c = 10^3     
    I_c₀ = 0
    S₀_c = N_c-I_c₀ 

    σₘ = 0.0158#0.0158

    #Calculate mutation matrix
    m,c_m = create_mutation_matrix(Nₚ,σₘ)
    
    S = [S₀_a,S₀_b,S₀_c]     
    I = [Infected([0,0,0],z_fun(i,Nₚ) ) for i in 1:Nₚ]  
    I[z_start_idx].number[1] = I_a₀
    I[z_start_idx].number[2] = I_b₀
    I[z_start_idx].number[3] = I_c₀
    @show I[z_start_idx]

    @show I_a_eq =  (1-c/(γ_a(μ_a)*K_aa))*N_a
    
    @show γ_a.(z_start)*K_aa/c
    @show γ_b.(z_start)*K_bb/c
    @show γ_c.(z_start)*K_cc/c
    
    #Plot R₀
    plot(z_vec,γ_a.(z_vec)*K_aa/c,xlabel="Strategy",ylabel="R₀",label="Species A") 
    plot!(z_vec,γ_b.(z_vec)*K_bb/c,label="Species B") 
    plot!(z_vec,γ_c.(z_vec)*K_cc/c,label="Species C")
    
    t₀ = 0.0
    t_end = 2000.0

    K = [K_aa K_ab K_ac;K_ab K_bb K_bc;K_ac K_bc K_cc]
    γ = (γ_a,γ_b,γ_c)
    N = [N_a,N_b,N_c]

    p = c,μₘ,σₘ,γ,K,N,c_m,t_end

    n_samples = 100

    S_vec,I_vec,t_vec = run_Gillespie!(S,I,t₀,p,n_samples)     
    
    
    draw_strain_distribution(t_vec,I_vec,"strain_evo")
    draw_strain_evolution_plot(t_vec,I_vec,t₀,t_end,"trait_vs_time")
    draw_compartments(t_vec,I_vec,S_vec,"compartments")    
end

run_ex9()