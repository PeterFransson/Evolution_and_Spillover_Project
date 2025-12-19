#stochastic individual-based two species evolutionary SIS-model simulation using 
#A hybrid τ-leap algorithm (Zhu et al. 2011)
#continuous pathogen strategy space

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
    ϵ::Real #0.01 ∼ 0.1 Zhu et al. 2011
    N_c::Integer #Threshold for critical number of individuals 
end

function Base.:+(a::Infected,b::Infected)
    length(a.number)==length(b.number)||error("length(a.number)!=length(b.number)")
    a.strategy≈b.strategy||println("Warning: a.strategy!≈b.strategy")
    return Infected(a.number+b.number,a.strategy)
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
τ_fun(r,rand_num) = -log(rand_num)/r
function inv_IDX(IDX,limit)
    temp = mod(IDX,limit)
    
    if temp==0
        return limit
    else
        return temp
    end
end
function get_event_species_strain(idx,n_events,n_species)
    IDX = idx
    event = inv_IDX(IDX,n_events)
    IDX = div(IDX-state,n_events)+1
    species = inv_IDX(IDX,n_species)
    pathogen_strains = div(IDX-species,n_species)+1

    return (event,species,pathogen_strains)
end

function update_state!(S::Vector{U},
    I::Vector{Infected},
    pathogen_strains::Integer,
    species::Integer,        
    event::Integer,
    number_of_events::Integer,
    σₘ::Real) where {U<:Integer,T<:Integer}

    if event==1 #recover
        I_z = I[pathogen_strains].number
        if sum(I_z)-1>0
            I[pathogen_strains].number[species] -= number_of_events
        else
            popat!(I,pathogen_strains)
        end
        S[species]+=number_of_events        
    elseif event==2 #mutation (always critical event)    
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
    elseif event==3 #infection
        I[pathogen_strains].number[species] += number_of_events
        S[species]-=number_of_events       
    else
        error("event error: sate either i, m, or r")
    end

    return nothing
end

function run_non_crit_simulations!(S::Vector{U},
    I::Vector{Infected},
    σₘ::Real,
    r_num,
    e_Σσ,
    n_e::Integer,
    τ::Real,
    n_events::Integer,
    n_species::Integer) where {U<:Integer}

    for i in 1:n_e 
        if e_Σσ[i]
            λ = r_num[i]*τ
            k = pois_rand(λ)
            if k>0
                event,species,pathogen_strains = get_event_species_strain(i,n_events,n_species) 
                update_state!(S,I,pathogen_strains,species,event,k,σₘ)    
            end       
        end
    end

    return nothing
end

function find_critical_event(r_num,e_Σσ,r_sum,rn_num)
    threshold = rn_num*r_sum  
    IDX = 0 
    temp_sum = 0   
    while temp_sum<threshold 
        IDX += 1 
        temp_sum += r_num[IDX]*!e_Σσ[IDX]     
    end 
    return IDX
end

function SI_Gillespie!(S::Vector{U},
    I::Vector{Infected},
    n_species::T,
    t,
    p) where {U<:Integer,T<:Integer}

    c,μₘ,σₘ,γ,K,N,N_crit,ϵ,t_end = p

    if isempty(I)==false
        n_pathogen_strains = length(I)   
        n_events = 3  
        
        get_IDX(k,j,i) = k+n_events*(j-1)+n_species*n_events*(i-1)
        
        n_e = n_pathogen_strains*n_species*n_events #Numer of events

        r_num = zeros(n_e) #Rates (propensities)
        e_Σσ = zeros(Bool,n_e) #partition events into critical (false) and non critical (true)
        r_sum = 0 #Sum rate of critical events        
        τ_min = Inf #Tau leap
        
        #Calculate rates
        for i in 1:n_pathogen_strains #Pathogen varriant loop
            I_z = I[i].number
            strategy_z = I[i].strategy
            for j in 1:n_species #Species loop                       
                 
                rᵣ = I_z[j]*c #k = 1 recovery
                rₘ = I_z[j]*μₘ #k = 2 virus mutation 
                rᵢ = S[j]*γ[j](strategy_z)*sum(K[j,:].*I_z./N) #k = 3 infection

                r_num[get_IDX(1,j,i)] = rᵣ
                r_num[get_IDX(2,j,i)] = rₘ
                r_num[get_IDX(3,j,i)] = rᵢ   
                
                r_sum += rₘ #Mutation events are always critical events 

                if I_z[j]>N_crit #Non-crit

                    e_Σσ[get_IDX(1,j,i)] = true #Non critical event
                    e_Σσ[get_IDX(3,j,i)] = true #Non critical event 

                    τ_temp = ϵ*I_z[j]/(rᵢ-rᵣ-rₘ) #τ estimate
                    if τ_temp<τ_min
                        τ_min=τ_temp                        
                    end   
                else
                    r_sum += rᵣ+rᵢ   
                end
            end       
        end

        Δt  = τ_fun(r_sum,rand()) #waiting time for critical event

        if Δt<τ_min #simulate a critical event
            rn_num = rand()
            IDX = find_critical_event(r_num,e_Σσ,r_sum,rn_num)            
            event,species,pathogen_strains = get_event_species_strain(IDX,n_events,n_species)
            update_state!(S,I,pathogen_strains,species,event,1,σₘ)
        end       

        h = min(Δt,τ_min)
        #simulate non critical events
        run_non_crit_simulations!(S,
        I,        
        σₘ,
        r_num,
        e_Σσ,
        n_e,
        h,
        n_events,
        n_species)  
        
        return t+h
    else
        println("Extinction")
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
            push!(S_vec,copy(S))         
            push!(t_vec,t)            
            if  length(I)>0                              
                push!(I_vec,deepcopy(I))                  
            else                
                push!(I_vec,Infected[])
            end
            sample_nr += 1            
        end       
    end 
        
    return (S_vec,I_vec,t_vec)
end

function get_evo_img_idx(strain::Real,strain_start::Real,strain_end::Real,Δstrain::Real)
    strain>strain_start||error("strain<=strain_start")
    strain<strain_end||error("strain>=strain_start")
    return floor(Integer,(strain-strain_start)/Δstrain)+1 #idx
end

function draw_evolution(t_vec,I_vec,img_name::String,z_start::Real,z_end::Real,nₚ::Integer)     
    n_time = length(t_vec)

    Δz = (z_end-z_start)/(nₚ-1)
    strategy_strain = [z_start+Δz*i for i in 0:nₚ-1]    

    ev_pl = zeros(n_time,nₚ) #evolution plot     

    for t_idx in 1:n_time
        if !isempty(I_vec[t_idx])
            I_strain = [sum(strain.number) for strain in I_vec[t_idx]]
            idxs = [get_evo_img_idx(strain.strategy,z_start,z_end,Δz) for strain in I_vec[t_idx]]
            
            I_tot = sum(I_strain)

            if I_tot>0                           
                ev_pl[t_idx,idxs] .= 1.0           
            end 
        end       
    end    
   
    plt = heatmap(strategy_strain,
    t_vec,
    ev_pl,
    c=cgrad(:grayC),
    cbar=false,    
    xlabel="strategy (Z)", ylabel="Time")
    
    savefig(plt,img_name*".svg")
end

function draw_compartments(t_vec,I_vec,S_vec,img_name::String)    
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
   
    savefig(plt,img_name*".svg")
end

function run_sample!(S,I,t₀,p,n_samples,file_name::AbstractString,syspar::SystemParameters,simpar::StocSimPar)     
    S_vec,I_vec,t_vec = run_Gillespie!(S,I,t₀,p,n_samples) 

    JLD2.@save file_name*".jld2" S_vec I_vec t_vec syspar simpar

    return nothing
end

function run_sample(file_name::AbstractString,syspar::SystemParameters,simpar::StocSimPar)    
    γ = syspar.γ_a #Recovery rate, assumption syspar.γ_a=syspar.γ_b   
    σ²,amplitude = syspar.σ²,syspar.τ_max
    μ_a = syspar.z_a
    μ_b = syspar.z_b
    
    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)    
   
    #Parameters
    #Maximum intraspecific transmission rate
    β_aa_max = syspar.c_aa
    β_bb_max = syspar.c_bb

    β_ab_max = syspar.c_ab
    
    μₘ = simpar.μₘ

    z_start = syspar.z_a

    #Initial states
    N_a = syspar.N_a     
    N_a>2||error("N_a<3")
    I_a₀ = ceil(Int,N_a*(1-γ/β_aa_max)) 
    S₀_a = N_a-I_a₀ 
    N_b = syspar.N_b  
    N_b>2||error("N_b<3") 
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 
    σₘ = simpar.σₘ
    
    S = [S₀_a,S₀_b]     
    I = [Infected([I_a₀,I_b₀],z_start)] 
        
    t₀ = simpar.t₀
    t_end = simpar.t_end

    β_max = [β_aa_max β_ab_max;β_ab_max β_bb_max]
    τ = (τ_a,τ_b)
    N = [N_a,N_b]

    p = γ,μₘ,σₘ,τ,β_max,N,simpar.N_c,simpar.ϵ,t_end
    
    n_samples = simpar.n_samples   

    run_sample!(S,I,t₀,p,n_samples,file_name,syspar,simpar)

    return nothing
end

function create_samples(folder_name::AbstractString,syspar::SystemParameters,simpar::StocSimPar;n_traj::Integer=8) 
    file_name = "sample"
          
    isdir(folder_name)||mkdir(folder_name)

    Threads.@threads for i in 1:n_traj 
        #Run simulation
        run_sample(folder_name*file_name*"_$(i)",syspar,simpar)
        println("Done: sample $(i)/$(n_traj)")
    end
    return nothing
end
function create_sample_fig(data_folder::AbstractString,figure_folder::AbstractString,n_traj::Integer,z_start::Real,z_end::Real,nₚ::Integer)
    file_name = "sample"  

    isdir(figure_folder)||mkdir(figure_folder)

    for i in 1:n_traj
        S_vec = JLD2.load(data_folder*file_name*"_$(i).jld2","S_vec")
        I_vec = JLD2.load(data_folder*file_name*"_$(i).jld2","I_vec")
        t_vec = JLD2.load(data_folder*file_name*"_$(i).jld2","t_vec")        

        figure_folder_temp = figure_folder*"/sample_$(i)/"
        isdir(figure_folder_temp)||mkdir(figure_folder_temp)

        draw_evolution(t_vec,I_vec,figure_folder_temp *"evol_plot",z_start,z_end,nₚ) 
        draw_compartments(t_vec,I_vec,S_vec,figure_folder_temp*"comp_plot")       
    end
end

function work_list()
    #--Create system parameter-- 

    #Basic parameters
    R₀_aa_max = 2.0
    R₀_bb_max = 2.0
    Δᵣ = 1.2
    c_ratio = 0.6

    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2    
    @show μ_b = μ_a+sqrt(2*σ²)*Δᵣ   
    
    z_min_lim = μ_a-0.05
    @show z_max_lim = μ_b+0.05
    
    c_crit = min(R₀_aa_max,R₀_bb_max)*2/(R₀_aa_max+R₀_bb_max)
    c = c_crit*c_ratio 
        
    N_a = 10^3
    N_b = 10^3    
            
    #Maximum intraspecific transmission rate
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
    β_ab_max = (β_aa_max+β_bb_max)/2*c

    syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude)

    #--Create stochastic simulation parameter-- 
    Nₚ = 300 #Number of strains
    μₘ = 0.03 #Mutation rate <----
    σₘ = 0.001#0.01 #0.0158 <----
    t₀ = 0.0 #<----
    t_end = 2500.0 #<----
    n_samples = 2500 #<----   
    ϵ = 0.04
    N_c = 10

    simpar = StocSimPar(Nₚ,μₘ,σₘ,t₀,t_end,n_samples,ϵ,N_c)  

    n_traj=5

    sub_folder_name = "branching_2_v2/"
    data_folder_name = "./output/AD_2_Species/stochastic_simulation/continuous_model/"*sub_folder_name
    figure_folder = "./fig/AD_2_Species/stochastic_simulation/continuous_model/"*sub_folder_name

    isdir(figure_folder)||mkdir(figure_folder)
    
    draw_PIP_n_coex_region(figure_folder,syspar)      

    #Create samples
    create_samples(data_folder_name,syspar,simpar;n_traj=n_traj)
    #Create sample figures
    create_sample_fig(data_folder_name,figure_folder,n_traj,z_min_lim,z_max_lim,2000)

    return nothing 
end

work_list()