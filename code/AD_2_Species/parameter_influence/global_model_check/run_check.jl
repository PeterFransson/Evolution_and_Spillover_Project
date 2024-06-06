#Single-strain-system R₀
function R₀_fun(Δz,Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)
    a  = exp(-Δz^2/(2*σ²))
    b =  exp(-(Δx^2-2*Δz*Δx)/(2*σ²))
    R₀_mean = (R₀_aa_max+R₀_bb_max*b)/2
    R₀_diff = (R₀_aa_max-R₀_bb_max*b)/2

    return a*(R₀_mean+sqrt(R₀_diff^2+R₀_ab_max^2*b))
end 

function R₀_fun(Δz,syspar::SystemParameters)
    c_aa,c_bb,c_ab = syspar.c_aa,syspar.c_bb,syspar.c_bb
    γ_a,γ_b = syspar.γ_a,syspar.γ_b   
    z_a,z_b = syspar.z_a,syspar.z_b 
    σ²,τ_max = syspar.σ²,syspar.τ_max 

    Δx = abs(z_b-z_a)
    
    R₀_aa_max = τ_max*c_aa/γ_a
    R₀_bb_max = τ_max*c_bb/γ_b
    R₀_ab_max = τ_max*c_ab/γ_a #Assuming γ_a=γ_b

    return R₀_fun(Δz,Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)
end

function find_min_system_R₀(syspar::SystemParameters;n_points::Integer=100)
    f(x) = R₀_fun(x,syspar) 
    
    z_a,z_b = syspar.z_a,syspar.z_b
    Δx = abs(z_b-z_a)

    Δz_rel = collect(range(0.0,stop=1.0,length=n_points))
    Δz = Δz_rel*Δx

    return minimum(f.(Δz))
end


function compute()    
    δz_vec = collect(range(0.05,stop=0.15,length=100)) #Distance between species in resource space 
    c_vec = [0.1 0.15 0.2 0.3 0.4 0.5] #Intraspecific transmission rate coefficient

    JLD2.@save "./output/AD_2_Species/parameter_influence/model_check/c_dz_info.jld2" c_vec δz_vec
    
    #maximum intraspecific basic reproduction number
    R₀_aa_max = 2.0 
    R₀_bb_max = 1.5 

    μ_a = 0.2 #Position of species A in resource space 
    γ = 0.1 #Recovery rate. We assume recovery rate (γ) is equal among species
    z_start = μ_a-0.02
    σ²,amplitude = 0.0025,1.0    
    Nₚ = 1500    

    #Parameters (Maximum intraspecific transmission rate)
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
        
    N_a = 10^3 #Population species A   
    N_b = N_a*1.0 #Populaation species B

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
        
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    t_end = 8000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]      
    Threads.@threads for i in eachindex(c_vec)
        for j in eachindex(δz_vec)
            δz =  δz_vec[j] #Distance between species in resource space 
            c = c_vec[i] #Intraspecific transmission rate coefficint

            file_name = "output_c_$(i)_dz_$(j)"
            
            β_ab_max = (β_aa_max+β_bb_max)/2*c #Maximum interspecific transmission rate
            μ_b = μ_a+δz #Position of species B in resource space         
            z_end = μ_b+0.02

            option = z_start,z_end,Nₚ,u₀,tspan

            syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude)

            strats = singular_strategies(syspar,option)

            min_sys_R0 = find_min_system_R₀(syspar)
            
            JLD2.@save "./output/AD_2_Species/parameter_influence/model_check/"*file_name*".jld2" strats syspar min_sys_R0
        end
    end    
end

function Base.:(==)(a::SingularStrat,b::SingularStrat)
    return a.conv_stable==b.conv_stable && a.evo_stable==b.evo_stable
end
function Base.:(==)(a::Vector{SingularStrat},b::Vector{SingularStrat})
    return length(a)==length(b) && prod(a .== b)==1    
end

function find_singular_strategies_classes()
    c_vec = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/c_dz_info.jld2","c_vec")
    δz_vec = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/c_dz_info.jld2","δz_vec")

    C_idx = 1

    strategy_δz = Tuple{Vector{SingularStrat},Real}[]
    strat_init = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/output_c_$(C_idx)_dz_$(1).jld2","strats")
    min_sys_R0_init = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/output_c_$(C_idx)_dz_$(1).jld2","min_sys_R0")
    
    push!(strategy_δz,(strat_init,min_sys_R0_init))    

    for j = eachindex(δz_vec)
        file_name = "output_c_$(C_idx)_dz_$(j)" 
        strat_temp = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/"*file_name*".jld2","strats")
        min_sys_R0_temp = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/"*file_name*".jld2","min_sys_R0")

        if first(last(strategy_δz)) != strat_temp
            push!(strategy_δz,(strat_temp,min_sys_R0_temp))
        end
    end    

    for strat in strategy_δz
        println(strat)
    end
end

function draw_PIP_select()
    folder_name = "./fig/PIP_plots/misc"
    name = "test"
    isdir(folder_name)||mkdir(folder_name)
    fig_name = folder_name*"/"*name  

    #Load system parameters
    syspar_name = "output_c_$(1)_dz_$(50)" 
    syspar = JLD2.load("./output/AD_2_Species/parameter_influence/model_check/"*syspar_name*".jld2","syspar") 

    #Initial values    
    N_a = 10^3 #Population species A   
    N_b = N_a*1.0 #Populaation species B
    
    Nₚ = 3000 #Number of strains
    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
        
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    t_end = 8000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    z_start,z_end = syspar.z_a-0.02,syspar.z_b+0.02

    option = z_start,z_end,Nₚ,u₀,tspan 

    draw_PIP(fig_name,syspar,option)    
end

#compute() 
#find_singular_strategies_classes()
draw_PIP_select()