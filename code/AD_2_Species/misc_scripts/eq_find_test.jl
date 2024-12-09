#Test new method to determine equilibrium solution for the

function get_system_par(strain,syspar::SystemParameters)
    τ_a = τ_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_b = τ_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)

    c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b, = syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.γ_a,syspar.γ_b,syspar.N_a,syspar.N_b

    return (τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
end
function get_system_par(strain,S,syspar::SystemParameters)
    (τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b) = get_system_par(strain,syspar)

    S_a,S_b = S[1],S[2]

    return (S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
end

function get_system_par_alt(strain,syspar::SystemParameters)
    Δx₁z = (syspar.z_a-strain)/sqrt(2*syspar.σ²)
    Δx₂z = (syspar.z_b-strain)/sqrt(2*syspar.σ²)
    R₀₁ = syspar.τ_max*syspar.c_aa/syspar.γ_a
    R₀₂ = syspar.τ_max*syspar.c_bb/syspar.γ_b
    R₀₁₂ = syspar.τ_max*syspar.c_ab/syspar.γ_a #We assume syspar.γ_a=syspar.γ_b
    N₁ = syspar.N_a
    N₂ = syspar.N_b
    return (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂)
end

function system_matrix_det(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    
    det = (S_a*τ_a*c_aa/N_a-γ_a)*(S_b*τ_b*c_bb/N_b-γ_b)-(S_a*τ_a*c_ab/N_b)*(S_b*τ_b*c_ab/N_a)    
    
    return det
end
function system_matrix_det(strain,S,syspar::SystemParameters)

    (S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b) = get_system_par(strain,S,syspar)
    

    det = system_matrix_det(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return det
end

function calc_S(x,strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = get_system_par_alt(strain,syspar)

    S₁ = N₁/(1-x*R₀₁₂*exp(-Δx₁z^2)/N₂)
    S₂ = N₂-x*(S₁*R₀₁*exp(-Δx₁z^2)/N₁-1)
    return (S₁,S₂)
end
function system_matrix_det_fun(x,strain,syspar::SystemParameters)
    
    det = system_matrix_det(strain,calc_S(x,strain,syspar),syspar)
    return det
end


function get_x_lim_par(strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = get_system_par_alt(strain,syspar)
    temp = N₂*exp(Δx₁z^2)/R₀₁₂
    temp₂ = 1-R₀₁₂*exp(-Δx₁z^2)-R₀₁*exp(-Δx₁z^2)
    x_min = temp*(temp₂-sqrt(temp₂^2+4*R₀₁₂*exp(-Δx₁z^2)))/2
    x_max = temp*(1-R₀₁*exp(-Δx₁z^2))
    return (x_min,x_max)
end
function get_x_lims(strain,syspar::SystemParameters)
    (x_min,x_max) = get_x_lim_par(strain,syspar::SystemParameters)
    return (x_min,min(0,x_max))
end

function get_poly_coeff(strain,syspar::SystemParameters)
    (Δx₁z,Δx₂z,R₀₁,R₀₂,R₀₁₂,N₁,N₂) = get_system_par_alt(strain,syspar)
    
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

#Function to calculate the discriminant of a third-degree polynomial ax³+bx²+cx+d
function calc_discriminant(a,b,c,d)
    return b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d
end
function calc_discriminant_depressed_cubic(p,q)
    return 4*p^3+27*q^2
end

function get_depressed_cubic_coef(a,b,c,d)
    p = (3*a*c-b^2)/(3*a^2)
    q = (2*b^3-9*a*b*c+27*a^2*d)/(27*a^3)
    return (p,q)
end

function calc_depressed_cubic_num(p,q)
    return q^2/4+p^3/27
end

function discriminant_test(p,q)
    return calc_depressed_cubic_num(p,q)>=0 #implies that the discriminant of the equation is negative
end

#Cardon's forumla: 
#Finds the (only) real root of the the Depressed cubic t^3+pt+q=0 
#where discriminant_test(p,q) is true 
function cubic_sol_one_real(p,q)
    depressed_cubic_num = calc_depressed_cubic_num(p,q)

    u₁ = -q/2+sqrt(depressed_cubic_num)
    u₂ = -q/2-sqrt(depressed_cubic_num)

    return cbrt(u₁)+cbrt(u₂)
end
#Finds the three real distinct root of the the Depressed cubic t^3+pt+q=0 
#where discriminant_test(p,q) is false
function cubic_sol_three_real(p,q)
    t₁ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3)
    t₂ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3-2*π*1/3)
    t₃ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3-2*π*2/3)

    return (t₁,t₂,t₃)
end
#Finds all real roots of the depressed cubic t^3+pt+q=0
function cubic_sol(p,q)
    depressed_cubic_num = calc_depressed_cubic_num(p,q)
    if depressed_cubic_num≈0
        if p≈0
            return 0
        else
            t = -3*q/(2*p)
            return (3*q/p,t,t)            
        end    
    elseif depressed_cubic_num>0
        return cubic_sol_one_real(p,q)
    else
        return cubic_sol_three_real(p,q)
    end
end

#Finds all real solutions to ax³+bx²+cx+d=0 when (a,b,c,d) are real
function cubic_sol(a,b,c,d)
    (p,q) = get_depressed_cubic_coef(a,b,c,d)
    
    t_sols = cubic_sol(p,q)
    return t_sols.-b/(3*a)
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

#2-species-1-strain system
function find_eq_test()
    #--Create system parameter-- 

    R₀_aa_max = 1.4
    γ = 0.1 #Recovery rate    
    σ²,amplitude = 0.0025,1.0
    μ_a = 0.2
    N_a = 10^3
    N_b = 10^3  
    R₀_ratio = 1.5
    Δᵣ = 1.2
    c_ratio = 0.2
    strain_ratio = 0.5
        
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

    strain = syspar.z_a+(syspar.z_b-syspar.z_a)*strain_ratio

    @show check_system_matrix(strain,syspar)    

    println("x₂= $(syspar.z_b)")

    @show (x_min,x_max) = get_x_lim_par(strain,syspar)
    @show (x_min,x_max) = get_x_lims(strain,syspar)
    @show (a,b,c,d) = get_poly_coeff(strain,syspar)

    f(x) = a*x^3+b*x^2+c*x+d

    xVec = range(x_min,stop=x_max,length=200)
    detVec = [system_matrix_det_fun(x,strain,syspar) for x in xVec] 
    polyVec = f.(xVec)    
   
    plot(xVec,detVec)
    plot!(xVec,polyVec)      
       
    x_sol = bisection(x->system_matrix_det_fun(x,strain,syspar),x_min,x_max)
    @show x_sol

    x_sol = bisection(f,x_min,x_max)
    @show x_sol
    
    @show cubic_sol(a,b,c,d)

    #x_sol = cubic_sol(a,b,c,d)
    #@show x_sol 

    #=
    (S₁,S₂) = calc_S(x_sol,strain,syspar)

    println("S₁=$(S₁),S₂=$(S₂)")

    @show R = system_matrix(strain,(S₁,S₂),syspar)
    @show I = [syspar.N_a-S₁;syspar.N_a-S₂]
    @show R*I
    =#

    #Initial values 
      
    I_a₀ = 1  
    S₀_a = syspar.N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = syspar.N_b-I_b₀  

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]     

    @show eq_point = find_eq(u₀,tspan,(strain,syspar))
    println("S₁=$(eq_point.S_a),S₂=$(eq_point.S_b)")
    

    @show eq_point = find_eq(strain,syspar)
    println("S₁=$(eq_point.S_a),S₂=$(eq_point.S_b)")

    #@time find_eq(u₀,tspan,(strain,syspar))
    #@time find_eq(strain,syspar)    
end

find_eq_test()