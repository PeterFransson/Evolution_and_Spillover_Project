#Find and categories singular strategies for PIP parameters

#-Code to find equilibrium points for the two-species-one-strain system

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

#Get limits for feasable solutions in terms of the system cubic polynomial variable x
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

    length(x_feasible)==1||error("More then one or no feasable solution exist")

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

    tr = R₀₁*Δx₁z+R₀₂*Δx₂z
    det = (R₀₁*R₀₂-R₀₁₂*R₀₁₂)*Δx₁z*Δx₂z

    return (tr+sqrt(tr^2-4*det))/2
end

function f_fun(strain,syspar::SystemParameters)
    eq_point = find_eq(strain,syspar)
    return calculate_select_grad(strain,(eq_point.S_a,eq_point.S_b),syspar) 
end

function PIP_gen_test()
    #--Create system parameter--     
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2

    N_a = 10^3
    N_b = 10^3    
    R₀_aa_max = 1.4        
    R₀_ratio = 1.0 #R₀_ratio∈(0,∞) 
    Δᵣ = 2.0 #Δᵣ∈(0,∞) Distance between species (Δᵣ 1 = sqrt(2)*σ)
    c_ratio = 0.2 #c_ratio∈(0,1) IMPORTANT! This should not be to low otherwise its equal to a disconnected system 
    
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

    #strain = syspar.z_a+(syspar.z_b-syspar.z_a)*strain_ratio

    strain_x_vec = range(0.0,stop=1.0,length=200)

    f(x) =  f_fun(syspar.z_a+(syspar.z_b-syspar.z_a)*x,syspar)
    R₀(x) = get_system_R₀(syspar.z_a+(syspar.z_b-syspar.z_a)*x,syspar)

    #=
    @show singular_strat = find_zeros(f, 0.0,1.0)

    if length(singular_strat)>0
        for x_strain in singular_strat
            strain = syspar.z_a+(syspar.z_b-syspar.z_a)*x_strain
            @show check_convergence(strain,syspar)
            @show find_sec_inv_fit(strain,syspar)  
        end
    end

    plot(strain_x_vec,f.(strain_x_vec))
    plot!([0.0,1.0],[0.0,0.0])=#

    plot(strain_x_vec,R₀.(strain_x_vec))
    plot!([0.0,1.0],[1.0,1.0])
end

PIP_gen_test()