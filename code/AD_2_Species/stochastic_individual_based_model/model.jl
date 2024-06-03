#stochastic individual-based two species evolutionary SIS-model simulation using Gillespie algorithm
#with simple within-host competition (discretized pathogen strategy space) [updated algorithm] 

mutable struct Infected
    number::Vector{Integer}  #Number of infected
    strategy::Real #Pathogen strategy
end

#Stochastic simulation parameter
mutable struct StocSimPar
    Nₚ::Integer #Number of strains
    μₘ::Real #Mutation rate
    σₘ::Real #Mutation standard deviation
    t₀::Real #Simualtion start time 
    t_end::Real #Simualtion end time 
    n_samples::Integer #Number of time samples
end

function Base.:+(a::Infected,b::Infected)
    length(a.number)==length(b.number)||error("length(a.number)!=length(b.number)")
    a.strategy≈b.strategy||println("Warning: a.strategy!≈b.strategy")
    return Infected(a.number+b.number,a.strategy)
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

        #if τ[species](mutant_strategy)>τ[species](resident_strategy) #No within-host competition
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

            #ev_pl[t_idx,:] .= 1.0.-q_strain
            ev_pl[t_idx,q_strain.>0.01] .= 0.0
        end       
        
    end    
   
    plt = heatmap(strategy_strain,
    t_vec,
    ev_pl,
    c=cgrad(:grays),
    xlabel="Time", ylabel="strategy (Z)")

    #savefig(plt,"./fig/AD_2_Species/stochastic_simulation/discrete_model/"*img_name*".svg")
    savefig(plt,img_name*".svg")
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
    #savefig(plt,"./fig/AD_2_Species/stochastic_simulation/discrete_model/"*img_name*".svg")
    savefig(plt,img_name*".svg")
end

function run_sample!(S,I,t₀,p,n_samples,file_name,syspar::SystemParameters,simpar::StocSimPar)     
    S_vec,I_vec,t_vec = run_Gillespie!(S,I,t₀,p,n_samples) 

    JLD2.@save "./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*".jld2" S_vec I_vec t_vec syspar simpar

    return nothing
end

function run_sample(file_name,syspar::SystemParameters,simpar::StocSimPar)    
    Nₚ = simpar.Nₚ #200 #Number of strains <---- 
    z_vec = z_fun.(1:Nₚ,Ref(Nₚ)) #Pathogen strains

    #maximum intraspecific basic reproduction number
    #R₀_aa_max = 2.44
    #R₀_bb_max = 2.44 

    γ = syspar.γ_a #assumption syspar.γ_a=syspar.γ_b #0.1 #Recovery rate    
    σ²,amplitude = syspar.σ²,syspar.τ_max#0.0025,1.0
    μ_a = syspar.z_a#0.2
    μ_b = syspar.z_b#0.33 
    #c = 0.8 #Intraspecific transmission rate coefficint

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)    
   
    #Parameters

    #Maximum intraspecific transmission rate
    @show β_aa_max = syspar.c_aa#R₀_aa_max*γ
    @show β_bb_max = syspar.c_bb#R₀_bb_max*γ

    @show β_ab_max = syspar.c_ab#(β_aa_max+β_bb_max)/2*c
    
    μₘ = simpar.μₘ#0.01 #Mutation rate <----

    @show z_start_idx = floor(Int64,μ_a*Nₚ)
    @show z_start = z_fun(z_start_idx,Nₚ) 

    #Initial states
    N_a = syspar.N_a#10^3     
    N_a>2||error("N_a<3")
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = syspar.N_b#10^3    
    N_b>2||error("N_b<3") 
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 
    σₘ = simpar.σₘ#0.0010 #0.0158 <----

    #Calculate mutation matrix
    m,c_m = create_mutation_matrix(Nₚ,σₘ)
    
    S = [S₀_a,S₀_b]     
    I = [Infected([0,0],z_fun(i,Nₚ) ) for i in 1:Nₚ]  
    I[z_start_idx].number[1] = I_a₀
    I[z_start_idx].number[2] = I_b₀
    @show I[z_start_idx]   
    
    t₀ = simpar.t₀#0.0 #<----
    t_end = simpar.t_end#250.0 #<----

    β_max = [β_aa_max β_ab_max;β_ab_max β_bb_max]
    τ = (τ_a,τ_b)
    N = [N_a,N_b]

    p = γ,μₘ,σₘ,τ,β_max,N,c_m,t_end

    n_samples = simpar.n_samples#1000 #<----   

    run_sample!(S,I,t₀,p,n_samples,file_name,syspar,simpar)

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

function infect_interpolation(t,I_vec,t_vec)
    t>t_vec[1]||return [sum(strain.number) for strain in I_vec[1]]
    t<t_vec[end]||return [sum(strain.number) for strain in I_vec[end]]
    t_end_idx = findfirst(t_vec.>t)
    t_start_idx = t_end_idx-1
    t₁ = t_vec[t_end_idx]
    t₀ = t_vec[t_start_idx]
    Δt = t₁-t₀
    a₀ = (t₁-t)/Δt
    a₁ = (t-t₀)/Δt
    return [sum(I_vec[t_start_idx][strain_idx].number)*a₀+sum(I_vec[t_end_idx][strain_idx].number)*a₁ for strain_idx in eachindex(I_vec[t_start_idx])]
end

function create_mean_traj(img_name,t_start,t_end,n_time,n_samples)

    file_name = "sample_" 

    I_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*"$(1).jld2","I_vec")
    n_strains = length(I_vec[1])
    strategy_strain = [strain.strategy for strain in I_vec[1]]

    t_vec_new = collect(range(t_start,stop=t_end,length=n_time)) #Time line for mean evolution plot    
    ev_pl = ones(n_time,n_strains) #evolution plot  
    
    #Calculate mean trajectory
    for sample in 1:n_samples 
        I_vec_temp = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*"$(sample).jld2","I_vec")
        t_vec_temp = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*file_name*"$(sample).jld2","t_vec")
        for (t_idx,t) in enumerate(t_vec_new)       
            ev_pl[t_idx,:]+=infect_interpolation(t,I_vec_temp,t_vec_temp)              
        end    
    end
    
    #Calcaulte frequency
    for t_idx in eachindex(t_vec_new)       
        
        I_tot = sum(ev_pl[t_idx,:])

        if I_tot>0
            q_strain = ev_pl[t_idx,:]/I_tot
            ev_pl[t_idx,:] .= 1.0.-q_strain
        end        
    end    
        
    plt = heatmap(strategy_strain,
    t_vec_new,
    ev_pl,
    c=cgrad(:grays),
    xlabel="Time", ylabel="strategy (Z)")

    savefig(plt,"./fig/AD_2_Species/stochastic_simulation/discrete_model/"*img_name*".svg")
end

function create_samples() 
    file_name = "sample"
    n_traj = 5 #Number of trajectories

    folder_name = "test_3/"  
    isdir("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name)||mkdir("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name)

    Threads.@threads for i in 1:n_traj

        #--Create system parameter-- 

        #maximum intraspecific basic reproduction number
        R₀_aa_max = 2.44
        R₀_bb_max = 2.44

        γ = 0.1 #Recovery rate    
        σ²,amplitude = 0.0025,1.0
        μ_a = 0.2
        μ_b = 0.33 
        c = 0.8
        
        N_a = 10^3
        N_b = 10^3    
                
        #Maximum intraspecific transmission rate
        β_aa_max = R₀_aa_max*γ
        β_bb_max = R₀_bb_max*γ
        β_ab_max = (β_aa_max+β_bb_max)/2*c

        syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude)

        #--Create stochastic simulation parameter-- 
        Nₚ = 300 #Number of strains
        μₘ = 0.01 #Mutation rate <----
        σₘ = 0.0001 #0.0158 <----
        t₀ = 0.0 #<----
        t_end = 1000.0 #<----
        n_samples = 1000 #<----   

        simpar=StocSimPar(Nₚ,μₘ,σₘ,t₀,t_end,n_samples)              

        #Run simulation
        run_sample(folder_name*file_name*"_$(i)",syspar,simpar)
        println("Done: sample $(i)/$(n_traj)")
    end
    return nothing
end
function create_sample_fig()
    #=
    t_start = 0.0
    t_end = 250.0
    n_time = 1000
    n_samples = 20
    create_mean_traj("mean_ev",t_start,t_end,n_time,n_samples)
    =#

    file_name = "sample"
    n_traj = 5 #Number of trajectories

    subfolder_name = "./fig/AD_2_Species/stochastic_simulation/discrete_model/"
    folder_name = "test_3/"  

    isdir(subfolder_name*folder_name)||mkdir(subfolder_name*folder_name)

    for i in 1:n_traj
        S_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name*file_name*"_$(i).jld2","S_vec")
        I_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name*file_name*"_$(i).jld2","I_vec")
        t_vec = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name*file_name*"_$(i).jld2","t_vec")
    
        syspar = JLD2.load("./output/AD_2_Species/stochastic_simulation/discrete_model/"*folder_name*file_name*"_$(i).jld2","syspar")       

        folder_name_temp = subfolder_name*folder_name*"/sample_$(i)/"
        isdir(folder_name_temp)||mkdir(folder_name_temp)

        draw_evolution(t_vec,I_vec,folder_name_temp*"evol_plot") 
        draw_compartments(t_vec,I_vec,S_vec,folder_name_temp*"comp_plot") 

        z_start = syspar.z_a-0.02
        z_start>0||error("z_start<0")
        z_end = syspar.z_b+0.02
        
        #Initial values
        Nₚ = 1000
        N_a = syspar.N_a    
        N_a>2||error("N_a<3")
        I_a₀ = 1  
        S₀_a = N_a-I_a₀ 
         
        N_b = syspar.N_b   
        N_b>2||error("N_b<3")
        I_b₀ = 1
        S₀_b = N_b-I_b₀  

        t_start = 0
        t_end = 7000
        tspan = (t_start,t_end)

        u₀ = [S₀_a I_a₀;S₀_b I_b₀]  

        option = z_start,z_end,Nₚ,u₀,tspan    
        draw_PIP(folder_name_temp*"PIP",syspar,option)            
    end
end


#run_model()
create_samples() 
create_sample_fig()