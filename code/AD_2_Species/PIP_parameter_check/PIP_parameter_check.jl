#Find and categories singular strategies for PIP parameters

function create_syspar(R₀_aa_max::Real,
    γ::Real,
    σ²::Real,
    amplitude::Real,
    μ_a::Real,
    N_a::Integer,
    N_b::Integer,
    R₀_ratio::Real,
    Δᵣ::Real,
    c_ratio::Real)     
    #--Create system parameter-- 

    #maximum intraspecific basic reproduction number    
    R₀_bb_max = R₀_aa_max*R₀_ratio

    μ_b = μ_a+sqrt(2*σ²)*Δᵣ
    
    c_crit = min(R₀_aa_max,R₀_bb_max)*2/(R₀_aa_max+R₀_bb_max)
    c = c_crit*c_ratio
        
      
            
    #Maximum intraspecific transmission rate
    β_aa_max = R₀_aa_max*γ
    β_bb_max = R₀_bb_max*γ
    β_ab_max = (β_aa_max+β_bb_max)/2*c

    syspar = SystemParameters(μ_a,μ_b,β_aa_max,β_bb_max,β_ab_max,N_a,N_b,γ,γ,σ²,amplitude) 

    return syspar 
end 

#Find equilibrium point and calculate second order derivative of invasion fitness at this point
function find_sec_inv_fit(strain::Real,syspar::SystemParameters)
    eq_point = find_eq(strain,syspar)    
    S = (eq_point.S_a,eq_point.S_b)
    
    ddλ_max = calculate_sec_inv_fit(strain, S, syspar)

    return ddλ_max
end
#Find equilibrium point and calculate selection gradiant at this point
function find_selectgrad(strain::Real,syspar::SystemParameters)
    eq_point = find_eq(strain,syspar)    
    S = (eq_point.S_a,eq_point.S_b)    
    
    dλ_max = calculate_select_grad(strain,S,syspar)

    return dλ_max
end

#Check if singular strategy strain is convergence stable
function check_convergence(strain::Real,syspar::SystemParameters;h::Real=10^-4)
        
    dλ_max_1 =  find_selectgrad(strain+h,syspar)
    dλ_max_2 =  find_selectgrad(strain-h,syspar)

    return (dλ_max_1-dλ_max_2)/(2*h)<0
end

function get_system_R₀(strain::Real,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = SysPar2PIPPar(strain,syspar)
    
    expΔx₁z = exp(-Δx₁z^2)
    expΔx₂z = exp(-Δx₂z^2)

    tr = R₀₁*expΔx₁z+R₀₂*expΔx₂z
    det = (R₀₁*R₀₂-R₀₁₂*R₀₁₂)*expΔx₁z*expΔx₂z

    return (tr+sqrt(tr^2-4*det))/2
end

#Calculate the selection gradiant Sᵣ´(m=r) 
function f_fun(strain,syspar::SystemParameters)
    eq_point = find_eq(strain,syspar)
    return calculate_select_grad(strain,(eq_point.S_a,eq_point.S_b),syspar) 
end

function classify_singular_strat(sing_strat::Real,syspar::SystemParameters;h::Real=10^-4,no_pts::Integer=100)
    conv_check = check_convergence(sing_strat,syspar;h=h)    
    evo_check = find_sec_inv_fit(sing_strat,syspar)

    global_check = false #True if evolutionarily stable point is a global fitness maximum
    if evo_check<0
        eq_point = find_eq(sing_strat,syspar) #Find endemic state of the resident strain
        f_temp(x) = calculate_inv_fit(syspar.z_a+(syspar.z_b-syspar.z_a)*x,(eq_point.S_a,eq_point.S_b),syspar) 
        x_sol = find_zeros(f_temp,0.0,1.0,no_pts=no_pts)

        if isempty(x_sol)||length(x_sol)==1
            global_check = true
        end
    end
    
    return SingularStrat(sing_strat,conv_check,evo_check<0,global_check)
end

function PIP_gen_n_classify()
    no_pts = 100 #Number points used to determine initial bracket intervals in find_zeros
    R₀_crit = 1.01
    n_samples = 10^2 #Number of Sobol samples in the parameter space
    
    #Setup SobolSeq
    dim = 6 #number of parameters 
    s = SobolSeq(dim)

    #--Create system parameter--     
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0 
    μ_a = 0.2 

    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    sub_folder_name = "Sobol_2025_05_12/"
    file_name = "sample"
    name = folder_name*sub_folder_name*"samples/"*file_name

    isdir(folder_name*sub_folder_name)||mkdir(folder_name*sub_folder_name)
    isdir(folder_name*sub_folder_name*"samples/")||mkdir(folder_name*sub_folder_name*"samples/") #Folder to save samples
        
    N = 1
    for i in 1:n_samples        
        x = next!(s)

        N_a = trunc(Int,10^2+x[1]*(10^4-10^2))
        N_b = trunc(Int,10^2+x[2]*(10^4-10^2))   
        R₀_aa_max = 1.0+x[3]*(4.0-1.0) 
        R₀_bb_max = 0.5+x[4]*(6.0-0.5)     
        R₀_ratio = R₀_bb_max/R₀_aa_max #R₀_ratio∈(0,∞) 
        Δᵣ = 0.1+x[5]*(3.0-0.1) #Δᵣ∈(0,∞) Distance between species (1 Δᵣ = sqrt(2)*σ)
        c_ratio = 0.01+x[6]*(0.8-0.01)#c_ratio∈(0,1) IMPORTANT! This should not be to low otherwise its equal to a disconnected system 
        
        syspar = create_syspar(R₀_aa_max,
        γ,
        σ²,
        amplitude,
        μ_a,
        N_a,
        N_b,
        R₀_ratio,
        Δᵣ,
        c_ratio)    

        x2strain(x) = syspar.z_a+(syspar.z_b-syspar.z_a)*x
        f(x) =  f_fun(x2strain(x),syspar)
        R₀(x) = get_system_R₀(x2strain(x),syspar)
        g(x) = R₀(x)-R₀_crit

        R₀_zeros_x = find_zeros(g, 0.0,1.0,no_pts=no_pts) #Find points where R₀=R₀_crit

        if isempty(R₀_zeros_x)||(length(R₀_zeros_x) == 2&&minimum(R₀_zeros_x)>0&&maximum(R₀_zeros_x)<1)
            if isempty(R₀_zeros_x) #System R₀>1 for all strains
                sing_strats_x = find_zeros(f, 0.0,1.0,no_pts=no_pts) #Find all evolutionary singular strategies 
                R₀_bool = true #True if system R₀>1 all strains
            else #System R₀<1 in a interval [minimum(R₀_zeros_x),maximum(R₀_zeros_x)]
                sing_strats_x = find_zeros(f, 0.0,minimum(R₀_zeros_x),no_pts=no_pts) #Find all evolutionary singular strategies  in the interval [0,minimum(R₀_zeros_x)]
                sing_strats_x_r = find_zeros(f, maximum(R₀_zeros_x),1.0,no_pts=no_pts) #Find all evolutionary singular strategies  in the interval [maximum(R₀_zeros_x),1.0]
                append!(sing_strats_x,sing_strats_x_r)
                R₀_bool = false #True if system R₀<1 for some strains
            end
            

            if isempty(sing_strats_x)==false
                sing_strats = x2strain.(sing_strats_x)
                strats = [classify_singular_strat(sing_strat,syspar) for sing_strat in sing_strats]    
                
                JLD2.@save name*"_$(N).jld2" strats R₀_bool syspar

                N += 1 
            end
        end 
    end   
end

function is_in_strats(col::Vector{Tuple{Vector{SingularStrat},Bool}},item::Tuple{Vector{SingularStrat},Bool})
    for item_temp in col
        if item_temp[1]==item[1]&&item_temp[2]==item[2]
            return true
        end
    end
    return false
end

function find_unique_strats()
    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    sub_folder_name = "Sobol_2025_05_12/"
    
    file_names = readdir(folder_name*sub_folder_name*"samples/")

    unique_strats = Tuple{Vector{SingularStrat},Bool}[]
    unique_strats_files = String[]
    
    for file_name in file_names
        name = folder_name*sub_folder_name*"samples/"*file_name

        strats = JLD2.load(name,"strats")
        R₀_bool = JLD2.load(name,"R₀_bool")

        if isempty(unique_strats)
            push!(unique_strats,(strats,R₀_bool))
            push!(unique_strats_files,file_name)
        else
            if is_in_strats(unique_strats,(strats,R₀_bool))==false
                push!(unique_strats,(strats,R₀_bool))
                push!(unique_strats_files,file_name)
            end
        end
    end

    println("Unique PIPs")
    for (i,strat) in enumerate(unique_strats)
        n_sing_strats = length(strat[1])
        println("PIP type $(i): N sing strat:$(n_sing_strats), R₀>1:$(strat[2])")
        println("--Number of singular strategies:$(n_sing_strats)")
        for strat in strat[1]
            println("---(ESS:$(strat.evo_stable),Conv:$(strat.conv_stable)),Global:$(strat.global_evo))")
        end
        println("File example: $(unique_strats_files[i])")
        println("-----")
    end

    #Save 
    temp_name = "list_of_unique/"
    isdir(folder_name*sub_folder_name*temp_name)||mkdir(folder_name*sub_folder_name*temp_name) #Folder to save list of unique samples
    JLD2.@save folder_name*sub_folder_name*temp_name*"list.jld2" unique_strats_files
end

function draw_uniques()
    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    sub_folder_name = "Sobol_2025_05_12/"
    temp_name = "list_of_unique/"

    img_folder =  "./fig/AD_2_Species/PIP_par_classify/"

    file_list = JLD2.load(folder_name*sub_folder_name*temp_name*"list.jld2","unique_strats_files")
    
    isdir(img_folder*sub_folder_name)||mkdir(img_folder*sub_folder_name)
    isdir(img_folder*sub_folder_name*"PIPs/")||mkdir(img_folder*sub_folder_name*"PIPs/") #Folder to save PIPs

    for file_name in file_list
        
        syspar_file_path = folder_name*sub_folder_name*"samples/"*file_name
        img_file_path = img_folder*sub_folder_name*"PIPs/"*file_name[1:end-5]
        
        draw_PIP(img_file_path::String,syspar_file_path::String)
    end
end

function test_global()
    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    sub_folder_name = "Sobol_2025_05_12/"
    temp_name = "list_of_unique/"

    img_folder =  "./fig/AD_2_Species/PIP_par_classify/"

    file_list = JLD2.load(folder_name*sub_folder_name*temp_name*"list.jld2","unique_strats_files")

    syspar_file_path = folder_name*sub_folder_name*"samples/"*file_list[1]
    println(syspar_file_path)
    syspar = JLD2.load(syspar_file_path,"syspar")
    strats = JLD2.load(syspar_file_path,"strats")
    sing_strat = strats[1].strategy
    x_sing_strat = (sing_strat-syspar.z_a)/(syspar.z_b-syspar.z_a)
    println("ESS = $(strats[1].evo_stable), Strat=$(sing_strat),X_strat=$(x_sing_strat)")
    
    @show eq_point = find_eq(sing_strat,syspar) #Find endemic state of the resident strain

    f_temp(x) = calculate_inv_fit(syspar.z_a+(syspar.z_b-syspar.z_a)*x,(eq_point.S_a,eq_point.S_b),syspar) 

    f_temp_alt(x) = maximum(eigvals(system_matrix(syspar.z_a+(syspar.z_b-syspar.z_a)*x,(eq_point.S_a,eq_point.S_b),syspar)))
    
    x = range(0.0,1.0,1000)
    y = f_temp.(x) 
    y_alt = f_temp_alt.(x)

    x_sol = find_zeros(f_temp,0.0,1.0,no_pts=100)
    @show x_sol

    plot(syspar.z_a.+(syspar.z_b-syspar.z_a)*x,y)
    plot!(syspar.z_a.+(syspar.z_b-syspar.z_a)*x,y_alt)
    #plot!([sing_strat,sing_strat],[0.0,0.1])

    #@show f_temp(x_sing_strat)
    #@show f_temp_alt(x_sing_strat)

    #x_sol = find_zeros(f_temp,0.0,1.0,no_pts=no_pts)
end

PIP_gen_n_classify()
find_unique_strats()
draw_uniques()
#test_global()