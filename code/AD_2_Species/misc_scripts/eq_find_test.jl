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

#2-species-1-strain system
function find_eq_test()
    #--Create system parameter-- 

    R₀_ratio = 1.5
    Δᵣ = 1.0
    c_ratio = 0.5
    strain = 0.21
    
    syspar = create_syspar(R₀_ratio,Δᵣ,c_ratio)
    println("x₂= $(syspar.z_b)")

    @show (x_min,x_max) = get_x_lim_par(strain,syspar)
    @show (x_min,x_max) = get_x_lims(strain,syspar)

    xVec = range(x_min,stop=x_max,length=200)
    detVec = [system_matrix_det_fun(x,strain,syspar) for x in xVec] 

    plot(xVec,detVec)

    x_sol = bisection(x->system_matrix_det_fun(x,strain,syspar),x_min,x_max)
    @show x_sol

    (S₁,S₂) = calc_S(x_sol,strain,syspar)

    @show R = system_matrix(strain,(S₁,S₂),syspar)
    @show I = [syspar.N_a-S₁;syspar.N_a-S₂]
    @show R*I
end

find_eq_test()