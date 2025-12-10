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
    #syspar.z_a <= strain <= syspar.z_b||error("$(syspar.z_a) <= $(strain) <= $(syspar.z_b)")    
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

#Find equilibrium point and calculate second order derivative of invasion fitness at this point
function find_sec_inv_fit(strain::Real,syspar::SystemParameters)
    eq_point = try_find_eq(strain,syspar)    
    S = (eq_point.S_a,eq_point.S_b)
    
    ddλ_max = calculate_sec_inv_fit(strain, S, syspar)

    return ddλ_max
end
#Find equilibrium point and calculate selection gradiant at this point
function find_selectgrad(strain::Real,syspar::SystemParameters)
    eq_point = try_find_eq(strain,syspar)    
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

#Finds the sub-intervals where the differentiable function f is positive (f>0.0) for
#a given closed interval
function find_pos_intervals(x_vals::Vector{R},f::Function,f_dif::Function;tol::Real=10^-4) where {R<:Real}
    #x_vals is ordered and includes the limit for a given intervals and 
    #the roots (f(x)=0) within in the interval
    x_intervals = Tuple{Real,Real}[] #Intervals where f>0
    n_x_vals = length(x_vals)

    n_x_vals>=2||error("length(x_vals)<2")

    isleft(x) = f(x)>0.0||(abs(f(x))<tol&&abs(f_dif(x))>tol&&f_dif(x)>0.0)
    isright(x) = f(x)>0.0||(abs(f(x))<tol&&abs(f_dif(x))>tol&&f_dif(x)<0.0)

    left_candidates = [isleft(x) for x in x_vals]
    right_candidates = [isright(x) for x in x_vals]
    
    left_idx = 0    
    
    for idx in 1:n_x_vals  
        if left_idx==0
            if left_candidates[idx]
                left_idx=idx
            end
        else
            if right_candidates[idx]
                push!(x_intervals,(x_vals[left_idx],x_vals[idx]))
                left_idx=0               
            end
        end        
    end

    return x_intervals
end

function find_pos_intervals(f::Function,f_dif::Function;tol::Real=10^-4,no_pts::Integer=100,digits::Integer=4)
    x_vals = [0.0,1.0]
    roots = find_zeros(f,0.0,1.0,no_pts=no_pts)
    append!(x_vals ,roots)
    unique!(x->round(x,digits=digits),x_vals)
    sort!(x_vals)    

    x_intervals = find_pos_intervals(x_vals,f,f_dif;tol=tol)
    return x_intervals
end

#Calculate the selection gradiant Sᵣ´(m=r) 
function f_fun(strain,syspar::SystemParameters)    
    eq_point = try_find_eq(strain,syspar)   
    return calculate_select_grad(strain,(eq_point.S_a,eq_point.S_b),syspar) 
end

function classify_singular_strat(sing_strat::Real,syspar::SystemParameters;h::Real=10^-4,no_pts::Integer=100,n_digits::Integer=4)
    x_sing_strat = (sing_strat-syspar.z_a)/(syspar.z_b-syspar.z_a)
    conv_check = check_convergence(sing_strat,syspar;h=h)    
    evo_check = find_sec_inv_fit(sing_strat,syspar)

    global_check = false #True if evolutionarily stable point is a global fitness maximum
    if evo_check<0        
        eq_point = try_find_eq(sing_strat,syspar) #Find endemic state of the resident strain
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

function get_singular_strats(syspar::SystemParameters;R₀_crit::Real=1.001,no_pts::Integer=100)
    x2strain(x) = syspar.z_a+(syspar.z_b-syspar.z_a)*x
    f(x) =  f_fun(x2strain(x),syspar)
    R₀(x) = get_system_R₀(x2strain(x),syspar)
    g(x) = R₀(x)-R₀_crit
    g_dif(x) = calculate_select_grad(x2strain(x),(syspar.N_a,syspar.N_b),syspar)

    R₀_zeros_intervals = find_pos_intervals(g,g_dif) #Find intervals where R₀>=R₀_crit

    !isempty(R₀_zeros_intervals)||println("Warning: no R₀>=R₀_crit found")
    length(R₀_zeros_intervals)<=2||println("Warning: more thant 2 regions where R₀>=R₀_crit found")

    if length(R₀_zeros_intervals)==1||length(R₀_zeros_intervals)==2
        if length(R₀_zeros_intervals)==1 #System R₀>1 for all strains
            l_val, r_val = R₀_zeros_intervals[1][1],R₀_zeros_intervals[1][2]

            sing_strats_x = find_zeros(f, l_val,r_val,no_pts=no_pts) #Find all evolutionary singular strategies 
            R₀_bool = true #True if system R₀>1 all strains
        else #System R₀<1 in a interval [r_val₁,l_val₂]
            l_val₁, r_val₁ = R₀_zeros_intervals[1][1],R₀_zeros_intervals[1][2]
            l_val₂, r_val₂ = R₀_zeros_intervals[2][1],R₀_zeros_intervals[2][2]

            sing_strats_x = find_zeros(f,l_val₁,r_val₁,no_pts=no_pts) #Find all evolutionary singular strategies  in the interval [0,minimum(R₀_zeros_x)]
            sing_strats_x_r = find_zeros(f,l_val₂,r_val₂,no_pts=no_pts) #Find all evolutionary singular strategies  in the interval [maximum(R₀_zeros_x),1.0]
            
            append!(sing_strats_x,sing_strats_x_r)
            R₀_bool = false #True if system R₀<1 for some strains
        end
        
        if isempty(sing_strats_x)==false
            sing_strats = x2strain.(sing_strats_x)
            strats = [classify_singular_strat(sing_strat,syspar) for sing_strat in sing_strats]               
        end
    end
    return (strats,R₀_bool)
end

function create_syspar(Δᵣ,c_ratio,R₀_aa_max,R₀_bb_max)
    #Create syspar
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0 #variance and amplitude τ-function
    μ_a = 0.2 #Position of species a in resource space     
    N_a = 1000
    N_b = 1000      
    R₀_ratio = R₀_bb_max/R₀_aa_max #R₀_ratio∈(0,∞) 
    #Δᵣ∈(0,∞) Distance between species (1 Δᵣ = sqrt(2)*σ)
    #c_ratio∈(0,1) IMPORTANT! This should not be to low otherwise its equal to a disconnected system 
    
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

    return syspar
end

function save_heatmap_samples(Δᵣ_vec::Vector{R},c_ratio_vec::Vector{R},R₀_aa_max::R,R₀_bb_max::R,filepath::String) where {R<:Real}
    n_Δᵣ = length(Δᵣ_vec)
    n_c_ratio = length(c_ratio_vec)

    for i in 1:n_c_ratio
        for j in 1:n_Δᵣ
            syspar = create_syspar(Δᵣ_vec[j],c_ratio_vec[i],R₀_aa_max,R₀_bb_max)
            strats,R₀_bool = get_singular_strats(syspar)  
            
            JLD2.@save filepath*"_$(i)_$(j).jld2" strats R₀_bool syspar
        end
    end

    return nothing
end

function calc_heatmap_mat(Δᵣ_vec::Vector{R},c_ratio_vec::Vector{R},filepath::String) where {R<:Real}
    generalist_strat_a = [SingularStrat(0.0,true,true,true)]
    generalist_strat_b = [SingularStrat(0.0,true,true,false)]

    evo_branching_strat = [SingularStrat(0.0,true,false,false)]

    bistable_strat_a = [SingularStrat(0.0,true,true,false),SingularStrat(0.0,false,false,false),SingularStrat(0.0,true,true,false)]
    bistable_strat_b = [SingularStrat(0.0,true,true,true),SingularStrat(0.0,false,false,false),SingularStrat(0.0,true,true,false)]
    bistable_strat_c = [SingularStrat(0.0,true,true,false),SingularStrat(0.0,false,false,false),SingularStrat(0.0,true,true,true)]
    bistable_strat_d = [SingularStrat(0.0,true,true,true),SingularStrat(0.0,false,false,false),SingularStrat(0.0,true,true,true)]

    R₀_strat_a = [SingularStrat(0.0,true,true,false),SingularStrat(0.0,true,true,false)]
    R₀_strat_b = [SingularStrat(0.0,true,true,true),SingularStrat(0.0,true,true,false)]
    R₀_strat_c = [SingularStrat(0.0,true,true,false),SingularStrat(0.0,true,true,true)]
    R₀_strat_d = [SingularStrat(0.0,true,true,true),SingularStrat(0.0,true,true,true)]

    n_Δᵣ = length(Δᵣ_vec)
    n_c_ratio = length(c_ratio_vec)

    heatmap_mat = zeros(Integer,(n_c_ratio,n_Δᵣ))

    for i in 1:n_c_ratio
        for j in 1:n_Δᵣ             
            
            strats = JLD2.load(filepath*"_$(i)_$(j).jld2","strats")
            R₀_bool = JLD2.load(filepath*"_$(i)_$(j).jld2","R₀_bool")

            if strats==generalist_strat_a||strats==generalist_strat_b
                heatmap_mat[i,j] = 1
            elseif strats==evo_branching_strat
                heatmap_mat[i,j] = 2
            elseif strats==bistable_strat_a||strats==bistable_strat_b||strats==bistable_strat_c||strats==bistable_strat_d
                heatmap_mat[i,j] = 3
            elseif R₀_bool==false&&(strats==R₀_strat_a||strats==R₀_strat_b||strats==R₀_strat_c||strats==R₀_strat_d)
                heatmap_mat[i,j] = 4    
            end            
        end        
    end

    return heatmap_mat
end

function draw_interspecific_influence_plot()
    Δᵣ_vec = collect(range(start=0.1,stop=2.5,length=200))
    c_ratio_vec = collect(range(start=0.01,stop=0.8,length=200))

    folders = ["paper_fig_R0a_2_0_R0b_1_1/",
    "paper_fig_R0a_2_0_R0b_1_2/",
    "paper_fig_R0a_2_0_R0b_1_5/",
    "paper_fig_R0a_2_0_R0b_1_7/",
    "paper_fig_R0a_2_0_R0b_2_0/",
    "paper_fig_R0a_2_0_R0b_2_5/"] 
    R₀_aa_maxs = [2.0,2.0,2.0,2.0,2.0,2.0]
    R₀_bb_maxs = [1.1,1.2,1.5,1.7,2.0,2.5]
        
    length(folders)==length(R₀_aa_maxs)==length(R₀_bb_maxs)||error("folders,R₀_aa_maxs,R₀_bb_maxs not the same size")
    n_jobs = length(folders)

    for i = 1:n_jobs
        R₀_aa_max = R₀_aa_maxs[i]
        R₀_bb_max = R₀_bb_maxs[i]
        
        foldername = "./output/AD_2_Species/interspecific_contact_influence/"*folders[i]
        imgfoldername = "./fig/AD_2_Species/interspecific_contact_influence/"*folders[i]

        filename = "sample"
        filepath = foldername*filename

        isdir(foldername)||(mkdir(foldername);save_heatmap_samples(Δᵣ_vec,c_ratio_vec,R₀_aa_max,R₀_bb_max,filepath))
        isdir(imgfoldername)||(mkdir(imgfoldername))

        
        heatmap_mat = calc_heatmap_mat(Δᵣ_vec,c_ratio_vec,filepath) 
        
        mygrad = cgrad([RGB(0.0,0.0,0.0),
        RGB(87/255,16/255,109/255),
        RGB(187/255,55/255,85/255),
        RGB(249/255,141/255,10/255),
        RGB(252/255,1.0,164/255)],         
        categorical = true)
        
        heatmap(Δᵣ_vec,c_ratio_vec,heatmap_mat,c=mygrad)        
        savefig(imgfoldername*"heatmap.svg")
        
        host_shift_heatmap_mat = calc_host_shift_heatmap_mat(Δᵣ_vec,c_ratio_vec,filepath)
        host_shift_gradient = cgrad([RGB(0.0,158/255,115/255),
        RGB(204/255,121/255,167/255),
        RGB(230/255,159/255,0.0)],
        categorical = true)
                 
        heatmap(Δᵣ_vec,c_ratio_vec,host_shift_heatmap_mat,color=host_shift_gradient)           
        savefig(imgfoldername*"host_shift_heatmap.svg")
    end
end

draw_interspecific_influence_plot()