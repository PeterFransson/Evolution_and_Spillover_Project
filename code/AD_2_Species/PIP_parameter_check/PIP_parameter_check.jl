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

function find_eq_numeric(strain,syspar::SystemParameters;no_pts::Integer=100)
    (a,b,c,d) = get_poly_coeff(strain,syspar)    
    (x_min,x_max) = get_x_lims(strain,syspar)

    f(x) = a*x^3+b*x^2+c*x+d

    x_feasible = find_zeros(f,x_min,x_max,no_pts=no_pts)

    length(x_feasible)==1||error("More then one or no feasible solution exists")

    (S₁,S₂) = calc_S(x_feasible[1],strain,syspar)

    eq_point = EqPoint(S₁,S₂,syspar.N_a-S₁,syspar.N_b-S₂)

    return eq_point
end

function try_find_eq(strain::Real,syspar::SystemParameters)   
    get_system_R₀(strain,syspar)>1.0||(return EqPoint(convert(typeof(0.0),syspar.N_a),convert(typeof(0.0),syspar.N_b),0.0,0.0)) 
    syspar.z_a <= strain <= syspar.z_b||error("$(syspar.z_a) <= $(strain) <= $(syspar.z_b)")    
    try         
        eq_point = find_eq(strain,syspar) #Analytical solution
        return eq_point 
    catch e_msg1  
        try
            eq_point = find_eq_numeric(strain,syspar) #Numerical solution
            return eq_point
        catch e_msg2
            JLD2.@save "./output/AD_2_Species/error/error_file.jld2" strain syspar
            error(e_msg2)            
        end             
    end    
end

#Calculate the selection gradiant Sᵣ´(m=r) 
function f_fun(strain,syspar::SystemParameters)
    #eq_point = find_eq(strain,syspar) 
    eq_point = try_find_eq(strain,syspar)   
    return calculate_select_grad(strain,(eq_point.S_a,eq_point.S_b),syspar) 
end

function classify_singular_strat(sing_strat::Real,syspar::SystemParameters;h::Real=10^-4,no_pts::Integer=100,n_digits::Integer=4)
    x_sing_strat = (sing_strat-syspar.z_a)/(syspar.z_b-syspar.z_a)
    conv_check = check_convergence(sing_strat,syspar;h=h)    
    evo_check = find_sec_inv_fit(sing_strat,syspar)

    global_check = false #True if evolutionarily stable point is a global fitness maximum
    if evo_check<0
        #eq_point = find_eq(sing_strat,syspar) #Find endemic state of the resident strain
        eq_point = try_find_eq(sing_strat,syspar)
        f_temp(x) = calculate_inv_fit(syspar.z_a+(syspar.z_b-syspar.z_a)*x,(eq_point.S_a,eq_point.S_b),syspar) 
        x_sol = find_zeros(f_temp,0.0,1.0,no_pts=no_pts)

        push!(x_sol,x_sing_strat)
        unique!(x->round(x,digits=n_digits),x_sol)

        if length(x_sol)==1
            global_check = true
        end
    end
    
    return SingularStrat(sing_strat,conv_check,evo_check<0,global_check)
end

function PIP_gen_n_classify(folder_name::String,sub_folder_name::String,options::Dict)
    no_pts = options["no_pts"] #Number points used to determine initial bracket intervals in find_zeros
    R₀_crit = options["R₀_crit"]
    n_samples = options["n_samples"] #Number of Sobol samples in the parameter space

    #folder_name = "./output/AD_2_Species/PIP_par_classify/"
    #sub_folder_name = "Sobol_2025_05_12/"

    R0_aa_min, R0_aa_max = options["R0_aa_min"],options["R0_aa_max"]
    R0_bb_min, R0_bb_max = options["R0_bb_min"],options["R0_bb_max"]
    N_a_min, N_a_max = options["N_a_min"],options["N_a_max"]
    N_b_min, N_b_max = options["N_b_min"],options["N_b_max"]
    Δᵣ_min, Δᵣ_max = options["Δᵣ_min"],options["Δᵣ_max"]
    c_ratio_min, c_ratio_max = options["c_ratio_min"],options["c_ratio_max"]
    
    #Setup SobolSeq
    dim = 6 #number of parameters 
    s = SobolSeq(dim)

    #--Create system parameter--     
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0 #variance and amplitude τ-function
    μ_a = 0.2 #Position of species a in resource space     
    
    file_name = "sample"
    name = folder_name*sub_folder_name*"samples/"*file_name

    isdir(folder_name*sub_folder_name)||mkdir(folder_name*sub_folder_name)
    isdir(folder_name*sub_folder_name*"samples/")||mkdir(folder_name*sub_folder_name*"samples/") #Folder to save samples
        
    N = 1
    for i in 1:n_samples        
        x = next!(s)

        N_a = trunc(Int,N_a_min+x[1]*(N_a_max-N_a_min))
        N_b = trunc(Int,N_b_min+x[2]*(N_b_max-N_b_min))   
        R₀_aa_max = R0_aa_min+x[3]*(R0_aa_max-R0_aa_min) 
        R₀_bb_max = R0_bb_min+x[4]*(R0_bb_max-R0_bb_min)     
        R₀_ratio = R₀_bb_max/R₀_aa_max #R₀_ratio∈(0,∞) 
        Δᵣ = Δᵣ_min+x[5]*(Δᵣ_max-Δᵣ_min) #Δᵣ∈(0,∞) Distance between species (1 Δᵣ = sqrt(2)*σ)
        c_ratio = c_ratio_min+x[6]*(c_ratio_max-c_ratio_min) #c_ratio∈(0,1) IMPORTANT! This should not be to low otherwise its equal to a disconnected system 
        
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
    for (idx,item_temp) in enumerate(col)
        if item_temp[1]==item[1]&&item_temp[2]==item[2]
            return (true,idx)
        end
    end
    return (false,nothing)
end

function find_unique_strats(folder_name::String,sub_folder_name::String)
    
    file_names = readdir(folder_name*sub_folder_name*"samples/")
    n_samples = length(file_names)

    unique_strats = Tuple{Vector{SingularStrat},Bool}[]
    unique_strats_files = String[]
    unique_strats_n = Integer[]
    
    for file_name in file_names
        name = folder_name*sub_folder_name*"samples/"*file_name

        strats = JLD2.load(name,"strats")
        R₀_bool = JLD2.load(name,"R₀_bool")

        if isempty(unique_strats)
            push!(unique_strats,(strats,R₀_bool))
            push!(unique_strats_files,file_name)
            push!(unique_strats_n,1)
        else
            unique_bool,unique_idx = is_in_strats(unique_strats,(strats,R₀_bool))
            if unique_bool==false
                push!(unique_strats,(strats,R₀_bool))
                push!(unique_strats_files,file_name)
                push!(unique_strats_n,1)
            else
                unique_strats_n[unique_idx] += 1 
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
        println("Number of occurrence: $(unique_strats_n[i])/$(n_samples)")
        println("-----")
    end

    #Save list of examples of unique PIPs
    temp_name = "list_of_unique/"
    isdir(folder_name*sub_folder_name*temp_name)||mkdir(folder_name*sub_folder_name*temp_name) #Folder to save list of unique samples
    JLD2.@save folder_name*sub_folder_name*temp_name*"list.jld2" unique_strats_files unique_strats_n
end

#Draw unique PIPs
function draw_uniques(folder_name::String,sub_folder_name::String,img_folder::String)
        
    temp_name = "list_of_unique/" #Folder pointing to list of examples of unique PIPs

    file_list = JLD2.load(folder_name*sub_folder_name*temp_name*"list.jld2","unique_strats_files")
    
    isdir(img_folder*sub_folder_name)||mkdir(img_folder*sub_folder_name)
    isdir(img_folder*sub_folder_name*"PIPs/")||mkdir(img_folder*sub_folder_name*"PIPs/") #Folder to save PIPs

    for file_name in file_list
        
        syspar_file_path = folder_name*sub_folder_name*"samples/"*file_name
        img_file_path = img_folder*sub_folder_name*"PIPs/"*file_name[1:end-5]
        
        draw_PIP(img_file_path,syspar_file_path)
    end
end

function run_work_list()
    #Folder names
    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    img_folder =  "./fig/AD_2_Species/PIP_par_classify/" 
    sub_folder_name = "Sobol_2025_05_27/" 
    
    isdir(folder_name*sub_folder_name)||mkdir(folder_name*sub_folder_name)
    isdir(folder_name*sub_folder_name*"options/")||mkdir(folder_name*sub_folder_name*"options/")
    
    options = Dict("no_pts" => 100, #Number points used to determine initial bracket intervals in the find_zeros algorithm to find singular strategies
    "R₀_crit" => 1.01,
    "n_samples" => 10^6, #Number of Sobol samples in the parameter space
    "R0_aa_min" => 1.0, 
    "R0_aa_max" => 6.0,
    "R0_bb_min" => 0.5, 
    "R0_bb_max" => 6.0, 
    "N_a_min" => 10^2, 
    "N_a_max" => 10^4,
    "N_b_min" => 10^2,
    "N_b_max" => 10^4,
    "Δᵣ_min" => 0.1, 
    "Δᵣ_max" => 3.5,
    "c_ratio_min" => 0.01,
    "c_ratio_max" => 0.8)

    JLD2.@save folder_name*sub_folder_name*"options/options.jld2" options

    PIP_gen_n_classify(folder_name,sub_folder_name,options)
    find_unique_strats(folder_name,sub_folder_name)
    draw_uniques(folder_name,sub_folder_name,img_folder)
end

function test_global()
    folder_name = "./output/AD_2_Species/PIP_par_classify/"
    sub_folder_name = "Sobol_2025_05_26/"
    temp_name = "list_of_unique/"

    img_folder =  "./fig/AD_2_Species/PIP_par_classify/"

    file_list = JLD2.load(folder_name*sub_folder_name*temp_name*"list.jld2","unique_strats_files")

    syspar_file_path = folder_name*sub_folder_name*"samples/"*file_list[3]
    println(syspar_file_path)
    syspar = JLD2.load(syspar_file_path,"syspar")
    strats = JLD2.load(syspar_file_path,"strats")
    sing_strat = strats[1].strategy
    x_sing_strat = (sing_strat-syspar.z_a)/(syspar.z_b-syspar.z_a)
    println("ESS = $(strats[1].evo_stable), Strat=$(sing_strat),X_strat=$(x_sing_strat)")
    
    @show eq_point = find_eq(sing_strat,syspar) #Find endemic state of the resident strain

    f_temp(x) = calculate_inv_fit(syspar.z_a+(syspar.z_b-syspar.z_a)*x,(eq_point.S_a,eq_point.S_b),syspar)     
    
    x = range(0.0,1.0,1000)
    y = f_temp.(x) 
    
    x_sol = find_zeros(f_temp,0.0,1.0,no_pts=100)
    @show x_sol    

    @show round(x_sol[1],digits=4)
    @show round(x_sing_strat,digits=4)   
    @show round(x_sol[1],digits=4)==round(x_sing_strat,digits=4)
    
    plot(syspar.z_a.+(syspar.z_b-syspar.z_a)*x,y)  
    
    #plot!([sing_strat,sing_strat],[0.0,0.1])

    #@show f_temp(x_sing_strat)
    #@show f_temp_alt(x_sing_strat)

    #x_sol = find_zeros(f_temp,0.0,1.0,no_pts=no_pts)
end

function test_error_file()
    file_name = "./output/AD_2_Species/error/error_file.jld2"

    @show syspar = JLD2.load(file_name ,"syspar")
    @show strain = JLD2.load(file_name,"strain")
    @show get_system_R₀(strain,syspar)

    @show Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂ = SysPar2PIPPar(strain,syspar)
    @show a,b,c,d = get_poly_coeff(strain,syspar)
    @show x_min,x_max = get_x_lims(strain,syspar)

    R₀_crit = 1.01
    x2strain(x) = syspar.z_a+(syspar.z_b-syspar.z_a)*x
    f(x) = a*x^3+b*x^2+c*x+d
    R₀(x) = get_system_R₀(x2strain(x),syspar)
    g(x) = R₀(x)-R₀_crit
    select_grad(x) =  f_fun(x2strain(x),syspar)
    
    @show a,b,c,d 
    @show (p,q) = get_depressed_cubic_coef(a,b,c,d)
    @show depressed_cubic_num = calc_depressed_cubic_num(p,q)
    #@show x_sols = cubic_sol(a,b,c,d)  
    #@show x_sol = bisection(f,x_min,x_max)
    @show N₁/R₀₁
    @show N₂/R₀₂    
    #@show S₁,S₂ = calc_S(x_sol,strain,syspar)
    
    xvec = range(x_min,stop=x_max,length=200)   
    norm_strain = range(0,stop=1.0,length=200) 
   
    yvec = f.(xvec)   
    plot(xvec,yvec) 

    plot(norm_strain,R₀.(norm_strain))

    #eq_point = find_eq(strain,syspar)
    eq_point =  try_find_eq(strain,syspar)

    @show eq_point.S_a 
    @show eq_point.S_b

    #@show R₀_zeros_x = find_zeros(g, 0.0,1.0,no_pts=100)
end 

#run_work_list()
#test_global()
test_error_file()