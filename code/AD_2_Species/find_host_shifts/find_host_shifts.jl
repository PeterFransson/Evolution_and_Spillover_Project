function find_host_shifts(folder_name::String,sub_folder_name::String)
    file_names = readdir(folder_name*sub_folder_name*"samples/")
    n_samples = length(file_names)

    n_Suitables = 0
    n_host_shifts = 0
    n_branching = 0 #Number of evolutionary branching
    n_generalist = 0 #Number of evolutionary generalist outcome

    for file_name in file_names
        name = folder_name*sub_folder_name*"samples/"*file_name

        strats = JLD2.load(name,"strats")
        syspar = JLD2.load(name,"syspar")

        R₀₁ = syspar.τ_max*syspar.c_aa/syspar.γ_a
        R₀₂ = syspar.τ_max*syspar.c_bb/syspar.γ_b   
        Δᵣx₂x₁ = (syspar.z_b-syspar.z_a)/sqrt(2*syspar.σ²)  
        R₀_x₁_in_2 = R₀₂*exp(-Δᵣx₂x₁^2) #R₀ in a pure host species 2 population
        
        if R₀₂>1.0&&R₀_x₁_in_2<1.0&&R₀₁>1.0

            n_Suitables += 1

            #Assume strats is ordered with repsect to strategy(position in resource space)            
            if  strats[1].evo_stable&&strats[1].conv_stable #ESS
                ESS_strain = strats[1].strategy
                Δx₂z = (ESS_strain-syspar.z_b)/sqrt(2*syspar.σ²)
                R₀_ESS_in_2 = R₀₂*exp(-Δx₂z^2)
                if R₀_ESS_in_2>1.0 
                    n_host_shifts += 1
                    if length(strats)==1 #Generalist outcome
                        n_generalist += 1
                    end
                end
            elseif !strats[1].evo_stable&&strats[1].conv_stable #Evolutionary Branching
                n_host_shifts += 1
                n_branching += 1
            end
        end        
    end

    println("Number of suitable paramter combinations: $(n_Suitables)")
    println("Number of host shift events: $(n_host_shifts)")
    println("Number of branching events: $(n_branching)")
    println("Number of generalist events: $(n_generalist)")

    result_df = DataFrame(n_Suitables = [n_Suitables],
    n_host_shifts = [n_host_shifts],
    n_branching = [n_branching],
    n_generalist = [n_generalist])
    CSV.write(folder_name*sub_folder_name*"host_shift_result.csv",result_df,delim=";")
    return nothing
end

function calc_host_shift_heatmap_mat(Δᵣ_vec::Vector{R},c_ratio_vec::Vector{R},filepath::String) where {R<:Real}

    n_Δᵣ = length(Δᵣ_vec)
    n_c_ratio = length(c_ratio_vec)

    heatmap_mat = zeros(Integer,(n_c_ratio,n_Δᵣ))

    for i in 1:n_c_ratio
        for j in 1:n_Δᵣ             
            
            strats = JLD2.load(filepath*"_$(i)_$(j).jld2","strats")
            syspar = JLD2.load(filepath*"_$(i)_$(j).jld2","syspar")

            R₀₁ = syspar.τ_max*syspar.c_aa/syspar.γ_a
            R₀₂ = syspar.τ_max*syspar.c_bb/syspar.γ_b   
            Δᵣx₂x₁ = (syspar.z_b-syspar.z_a)/sqrt(2*syspar.σ²)  
            R₀_x₁_in_2 = R₀₂*exp(-Δᵣx₂x₁^2) #R₀ in a pure host species 2 population
            
            if R₀₂>1.0&&R₀_x₁_in_2<1.0&&R₀₁>1.0            
                heatmap_mat[i,j] = 1
                #Assume strats is ordered with repsect to strategy(position in resource space)            
                if  strats[1].evo_stable&&strats[1].conv_stable #ESS
                    ESS_strain = strats[1].strategy
                    Δx₂z = (ESS_strain-syspar.z_b)/sqrt(2*syspar.σ²)
                    R₀_ESS_in_2 = R₀₂*exp(-Δx₂z^2)
                    if R₀_ESS_in_2>1.0 
                        heatmap_mat[i,j] = 2                
                    end
                elseif !strats[1].evo_stable&&strats[1].conv_stable #Evolutionary Branching
                    heatmap_mat[i,j] = 2
                end
            end               
        end        
    end

    return heatmap_mat
end