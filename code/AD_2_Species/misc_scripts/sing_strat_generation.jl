function create_syspar(R₀_ratio::Real,Δᵣ::Real,c_ratio::Real)     
    #--Create system parameter-- 

    #maximum intraspecific basic reproduction number
    R₀_aa_max = 1.4
    R₀_bb_max = R₀_aa_max*R₀_ratio

    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2
    μ_b = μ_a+sqrt(2*σ²)*Δᵣ
    
    c_crit = min(R₀_aa_max,R₀_bb_max)*2/(R₀_aa_max+R₀_bb_max)
    c = c_crit*c_ratio
        
    N_a = 10^3
    N_b = 10^3    
            
    #Maximum intraspecific transmission rate
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
    β_ab_max = (β_aa_max+β_bb_max)/2*c

    syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude) 

    return syspar 
end 

function strat_work_list()
    R₀_ratio_vec = [1.0,1.1,1.4,1.5,2.0]#[1.3]#[1.0,1.5,2.0]
    c_ratio_vec = [0.8,0.6,0.3,0.2,0.1]#[0.8]#[0.8,0.6,0.2]
    Δᵣ_vec = collect(range(1.0,stop=3.0,length=30))   
    
    for i in eachindex(R₀_ratio_vec)
        for j in eachindex(Δᵣ_vec)
            Threads.@threads for k in eachindex(c_ratio_vec)
                R₀_ratio = R₀_ratio_vec[i]               
                Δᵣ = Δᵣ_vec[j]                
                c_ratio = c_ratio_vec[k]                

                syspar = create_syspar(R₀_ratio,Δᵣ,c_ratio)  

                (strats,R₀_bool) = get_singular_strategies(syspar,syspar.z_a-0.02,syspar.z_b+0.02) #R₀_bool true if there is a region with R₀<1
                
                JLD2.@save "./output/AD_2_Species/parameter_influence/test/sample_R_$(i)_D_$(j)_C_$(k).jld2" strats R₀_bool syspar
            end            
        end        
    end

    return nothing
end

function isequal(a::Tuple{Vector{SingularStrat},Bool},b::Tuple{Vector{SingularStrat},Bool})
    return first(a)==first(b) && last(a)==last(b)    
end

function get_strat_types()
    R₀_ratio_vec = [1.0,1.1,1.4,1.5,2.0]#[1.3]#[1.0,1.5,2.0]
    c_ratio_vec = [0.8,0.6,0.3,0.2,0.1]#[0.8]#[0.8,0.6,0.2]
    Δᵣ_vec = collect(range(1.0,stop=3.0,length=30))   
    
    strat_types = Tuple{Vector{SingularStrat},Bool}[]
    
    for i in eachindex(R₀_ratio_vec)
        for j in eachindex(Δᵣ_vec)
            for k in eachindex(c_ratio_vec)
                strats = JLD2.load("./output/AD_2_Species/parameter_influence/test/sample_R_$(i)_D_$(j)_C_$(k).jld2","strats")
                R₀_bool = JLD2.load("./output/AD_2_Species/parameter_influence/test/sample_R_$(i)_D_$(j)_C_$(k).jld2","R₀_bool")
                
                if isempty(strat_types)
                    push!(strat_types,(strats,R₀_bool))
                end  
                if !any(isequal.(strat_types,Ref((strats,R₀_bool)))) 
                    push!(strat_types,(strats,R₀_bool))
                end         
            end
        end
    end

    for strat in strat_types
        println(strat)
    end
    
    return nothing
end

function check_single_strat_type()
    R₀_ratio_vec = [1.0,1.1,1.4,1.5,2.0]#[1.3]#[1.0,1.5,2.0]
    c_ratio_vec = [0.8,0.6,0.3,0.2,0.1]#[0.8]#[0.8,0.6,0.2]
    Δᵣ_vec = collect(range(1.0,stop=3.0,length=30))   

    target_strat = (SingularStrat[SingularStrat(NaN, true, false),
    SingularStrat(NaN, false, false),
    SingularStrat(NaN, true, false)],
    false)   
    
    for i in eachindex(R₀_ratio_vec)
        for j in eachindex(Δᵣ_vec)
            for k in eachindex(c_ratio_vec)
                strats = JLD2.load("./output/AD_2_Species/parameter_influence/test/sample_R_$(i)_D_$(j)_C_$(k).jld2","strats")
                R₀_bool = JLD2.load("./output/AD_2_Species/parameter_influence/test/sample_R_$(i)_D_$(j)_C_$(k).jld2","R₀_bool")
                
                if isequal(target_strat,(strats,R₀_bool))
                    println("i: $(i), j: $(j), k: $(k)")
                end                           
            end
        end
    end  
    
    return nothing
end

#strat_work_list()
#get_strat_types()
check_single_strat_type()