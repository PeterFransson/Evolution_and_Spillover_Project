#stochastic individual-based two species evolutionary SIS-model simulation using Gillespie algorithm
#with simple within-host competition (discretized pathogen strategy space) [updated algorithm] 

mutable struct Infected
    number::Vector{Integer}  #Number of infected
    strategy::Real #Pathogen strategy
end

τ_fun(r,rand_num) = -log(rand_num)/r

#Simpson's 1/3 rule
simpson13(f,a,b) = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6
#Discretize a continuous positive function f
function cont2disc(f,disc_samples)
    n = length(disc_samples)
    areas = simpson13.(Ref(f),disc_samples[1:n-1],disc_samples[2:n])
    area_tot = sum(areas)
    return areas/area_tot
end

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

function SIS_Gillespie!(S::Vector{U},
    I::Vector{Infected},
    n_species::T,
    t,
    p) where {U<:Integer,T<:Integer}

    γ,μₘ,σₘ,τ,K,N,c_m,t_end = p #Parameters
    
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
                
            rᵣ = I_z[j]*γ
            rₘ = I_z[j]*μₘ
            rᵢ = S[j]*τ[j](strategy_z)*sum(K[j,:].*I_z./N)

            τ_val = [τ_fun(rᵣ,r_num[i,j,1]),τ_fun(rₘ,r_num[i,j,2]),τ_fun(rᵢ,r_num[i,j,3])]
            state_temp = argmin(τ_val)
            if τ_val[state_temp]<Δt
                Δt = τ_val[state_temp]
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

        #if τ[species](mutant_strategy)>τ[species](resident_strategy)
            I[pathogen_strains_mutant].number[species]+=1
            I[pathogen_strains].number[species] -= 1                
        #end
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
        t = SIS_Gillespie!(S,I,n_species,t,p)    

        if t>t₀+sample_nr*Δt_sample
            #println("$(sample_nr)/$(n_samples)")
            push!(S_vec,copy(S))         
            push!(t_vec,t) 
            push!(I_vec,deepcopy(I))          
            sample_nr += 1            
        end       
    end 
        
    return (S_vec,I_vec,t_vec)
end

function draw_evolution(t_vec,I_vec,img_name)     
    n_time = length(t_vec)
    n_strains = length(I_vec[1])
    strategy_strain = [strain.strategy for strain in I_vec[1]]

    ev_pl = ones(n_time,n_strains) #evolution plot      

    for t_idx in 1:n_time

        I_strain = [sum(strain.number) for strain in I_vec[t_idx]]
        
        I_tot = sum(I_strain)

        if I_tot>0
            q_strain = I_strain/I_tot

            ev_pl[t_idx,:] .= 1.0.-q_strain
        end       
        
    end    
   
    plt = heatmap(strategy_strain,
    t_vec,
    ev_pl,
    c=cgrad(:grays),
    xlabel="Time", ylabel="strategy (Z)")

    savefig(plt,"./fig/AD_2_Species/stochastic_simulation/discrete_model/"*img_name*".svg")
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
    savefig(plt,"./fig/AD_2_Species/stochastic_simulation/discrete_model/"*img_name*".svg")
end

function run_sample!(S,I,t₀,p,n_samples,file_name)     
    S_vec,I_vec,t_vec = run_Gillespie!(S,I,t₀,p,n_samples) 

    JLD2.@save "./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*".jld2" S_vec I_vec t_vec

    return nothing
end

function run_sample(file_name)    
    Nₚ = 200 #Number of strains 
    z_vec = z_fun.(1:Nₚ,Ref(Nₚ)) #Pathogen strains

    #maximum intraspecific basic reproduction number
    R₀_aa_max = 2.44
    R₀_bb_max = 2.44 

    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2
    μ_b = 0.33 
    c = 0.8 #Intraspecific transmission rate coefficint

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)    
   
    #Parameters

    #Maximum intraspecific transmission rate
    @show β_aa_max = R₀_aa_max*γ
    @show β_bb_max = R₀_bb_max*γ

    @show β_ab_max = (β_aa_max+β_bb_max)/2*c
    
    μₘ = 0.001 #Mutation rate

    @show z_start_idx = floor(Int64,μ_a*Nₚ)
    @show z_start = z_fun(z_start_idx,Nₚ) 

    #Initial states
    N_a = 10^3     
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 
    σₘ = 0.0010 #0.0158

    #Calculate mutation matrix
    m,c_m = create_mutation_matrix(Nₚ,σₘ)
    
    S = [S₀_a,S₀_b]     
    I = [Infected([0,0],z_fun(i,Nₚ) ) for i in 1:Nₚ]  
    I[z_start_idx].number[1] = I_a₀
    I[z_start_idx].number[2] = I_b₀
    @show I[z_start_idx]   
    
    t₀ = 0.0
    t_end = 1000.0

    β_max = [β_aa_max β_ab_max;β_ab_max β_bb_max]
    τ = (τ_a,τ_b)
    N = [N_a,N_b]

    p = γ,μₘ,σₘ,τ,β_max,N,c_m,t_end

    n_samples = 1000

    #file_name = "run_1"

    run_sample!(S,I,t₀,p,n_samples,file_name)

    return nothing
end

function run_model()
    file_name = "run_1"

    run_sample(file_name) 

    S_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*".jld2","S_vec")
    I_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*".jld2","I_vec")
    t_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*".jld2","t_vec")

    draw_evolution(t_vec,I_vec,"evol_plot") 
    draw_compartments(t_vec,I_vec,S_vec,"comp_plot") 
end


function create_samples() 
    file_name = "sample"
    n_samples = 20
    Threads.@threads for i in 1:n_samples
        run_sample(file_name*"_$(i)")
        println("Done: sample $(i)/$(n_samples)")
    end
    return nothing
end

#run_model()
create_samples() 