#Find and categories singular strategies for PIP parameters

#-Start: Code to find equilibrium points for the two-species-one-strain system

function SysPar2PIPPar(strain,syspar::SystemParameters)
    Δx₁z = (syspar.z_a-strain)/sqrt(2*syspar.σ²)
    Δx₂z = (syspar.z_b-strain)/sqrt(2*syspar.σ²)
    R₀₁ = syspar.τ_max*syspar.c_aa/syspar.γ_a
    R₀₂ = syspar.τ_max*syspar.c_bb/syspar.γ_b
    R₀₁₂ = syspar.τ_max*syspar.c_ab/syspar.γ_a #We assume syspar.γ_a=syspar.γ_b
    N₁ = syspar.N_a
    N₂ = syspar.N_b
    return (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂)
end

#Get coefficients for the system cubic polynomial f(x) = ax³+bx²+cx+d. 
#The roots of f(x) are equilibrium points candidates
function get_poly_coeff(strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = SysPar2PIPPar(strain,syspar)
    
    expΔx₁z = exp(-Δx₁z^2)
    expΔx₂z = exp(-Δx₂z^2)    

    a_eq = (R₀₁*R₀₂-R₀₁₂*R₀₁₂)*expΔx₁z*expΔx₂z/(N₁*N₂)
    b_eq = -R₀₁*expΔx₁z/N₁
    c_eq = -R₀₂*expΔx₂z/N₂

    a_s1 = -N₁*R₀₁₂*expΔx₁z/N₂
    b_s1 = N₁

    a_1 = (R₀₁₂*expΔx₁z/N₂)^2 
    b_1 = -2*R₀₁₂*expΔx₁z/N₂
    c_1 = 1

    a_s1s2 = -N₁*R₀₁₂*expΔx₁z/N₂
    b_s1s2 = N₁*(1-R₀₁*expΔx₁z-R₀₁₂*expΔx₁z) 
    c_s1s2 = N₁*N₂

    a_s2 = (R₀₁₂*expΔx₁z/N₂)^2 
    b_s2 = -R₀₁₂*expΔx₁z*(2-R₀₁*expΔx₁z-R₀₁₂*expΔx₁z)/N₂
    c_s2 = (1-R₀₁*expΔx₁z-2*R₀₁₂*expΔx₁z)
    d_s2 = N₂

    a = (a_s2*c_eq)
    b = (a_s1s2*a_eq+b_s2*c_eq+a_1)    
    c = (b_s1s2*a_eq+a_s1*b_eq+c_s2*c_eq+b_1)
    d = (c_s1s2*a_eq+b_s1*b_eq+d_s2*c_eq+c_1)
    
    return (a,b,c,d)
end

#Get limits for feasible solutions in terms of the system cubic polynomial variable x
function get_x_lim_vals(strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = SysPar2PIPPar(strain,syspar)
    temp = N₂*exp(Δx₁z^2)/R₀₁₂
    temp₂ = 1-R₀₁₂*exp(-Δx₁z^2)-R₀₁*exp(-Δx₁z^2)
    x_min = temp*(temp₂-sqrt(temp₂^2+4*R₀₁₂*exp(-Δx₁z^2)))/2
    x_max = temp*(1-R₀₁*exp(-Δx₁z^2))
    return (x_min,x_max)
end
function get_x_lims(strain,syspar::SystemParameters)
    (x_min,x_max) = get_x_lim_vals(strain,syspar::SystemParameters)
    return (x_min,min(0,x_max))
end

#Calcualte the number of susceptibles (S₁,S₂) from the system cubic polynomial variable x
function calc_S(x,strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = SysPar2PIPPar(strain,syspar)

    S₁ = N₁/(1-x*R₀₁₂*exp(-Δx₁z^2)/N₂)
    S₂ = N₂-x*(S₁*R₀₁*exp(-Δx₁z^2)/N₁-1)
    return (S₁,S₂)
end

#Function to find the equilibrium point for the two-species-one-pathogen system when the system R₀>0 and all 
#parameters are >0
function find_eq(strain,syspar::SystemParameters)
    (a,b,c,d) = get_poly_coeff(strain,syspar)
    x_sols = cubic_sol(a,b,c,d)  
    (x_min,x_max) = get_x_lims(strain,syspar)

    x_feasible = []
    for x_sol in x_sols      
        if x_min<x_sol<x_max
            push!(x_feasible,x_sol)
        end
    end

    length(x_feasible)==1||error("More then one or no feasible solution exists")

    (S₁,S₂) = calc_S(x_feasible[1],strain,syspar)

    eq_point = EqPoint(S₁,S₂,syspar.N_a-S₁,syspar.N_b-S₂)

    return eq_point
end
#-End: Code to find equilibrium points for the two-species-one-strain system

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

        if isempty(x_sol)||length(x_sol)<2
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
    file_name = "Sample"
    name = folder_name*sub_folder_name*file_name

    isdir(folder_name*sub_folder_name)||mkdir(folder_name*sub_folder_name)
        
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
    
    file_names = readdir(folder_name*sub_folder_name)

    unique_strats = Tuple{Vector{SingularStrat},Bool}[]

    for file_name in file_names
        name = folder_name*sub_folder_name*file_name

        strats = JLD2.load(name,"strats")
        R₀_bool = JLD2.load(name,"R₀_bool")

        if isempty(unique_strats)
            push!(unique_strats,(strats,R₀_bool))
        else
            if is_in_strats(unique_strats,(strats,R₀_bool))==false
                push!(unique_strats,(strats,R₀_bool))
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
        println("-----")
    end
end

PIP_gen_n_classify()
find_unique_strats()