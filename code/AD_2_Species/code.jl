#Simple SIS-model example with two species

#Data type to store parameters for the two-species system
mutable struct SystemParameters
    z_a::Real #Position of species A in resource space [0,1] 
    z_b::Real #Position of species B in resource space [0,1] 
    c_aa::Real #Within population contact rate of species A
    c_bb::Real #Within population contact rate of species B
    c_ab::Real #Between population contact rate of species A and B
    N_a::Real #Population of species A
    N_b::Real #Population of species B
    γ_a::Real #Recovery rate of species A
    γ_b::Real #Recovery rate of species B
    σ²::Real #Niche width
    τ_max::Real #Maxium τ value (probability of successful infection)
end

#Save parameters CSV-file
function savesystemparameters(p::SystemParameters,filename::String)
    param_df = DataFrame(z_a = [p.z_a],
    z_b = [p.z_b],
    c_aa = [p.c_aa],
    c_bb = [p.c_bb],
    c_ab = [p.c_bb],
    N_a = [p.N_a],
    N_b = [p.N_b],
    gamma_a = [p.γ_a],
    gamma_b = [p.γ_b],
    sigma_2 = [p.σ²],
    tau_max = [p.τ_max])
    CSV.write(filename*".csv",param_df,delim=";")
    return nothing
end

#Data type stores equilibrium point for the two-species-one-pathogen system
mutable struct EqPoint{T<:Real}
    S_a::T
    S_b::T
    I_a::T
    I_b::T  
end

#Struct containing information of a singular strategy
struct SingularStrat
    strategy::Real #Position in compability space of evolutionarily singular strategy
    conv_stable::Bool #True if the strategy is convergence stable
    evo_stable::Bool #True if the strategy is evolutionarily stable
    global_evo::Bool #True if the strategy is evolutionarily stable and a global fitness maximum
end

function Base.:(==)(a::SingularStrat,b::SingularStrat)
    return a.conv_stable==b.conv_stable && a.evo_stable==b.evo_stable && a.global_evo==b.global_evo
end
function Base.:(==)(a::Vector{SingularStrat},b::Vector{SingularStrat})
    return length(a)==length(b) && prod(a .== b)==1    
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
τ_fun(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))
τ_fun(x,μ,σ²,amplitude) = amplitude*τ_fun(x,μ,σ²)
τ_prime_fun(x,μ,σ²,amplitude) = -amplitude*exp(-(x-μ)^2/(σ²*2))*2*(x-μ)/(σ²*2) #derivative of τ_fun
τ_d_prime_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))*4*(x-μ)^2/(σ²*2)^2-amplitude*exp(-(x-μ)^2/(σ²*2))*2/(σ²*2) #Second order derivative of τ_fun

#Calculate invasion fitness
function calculate_inv_fit(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    r_11 = S_a*τ_a*c_aa/N_a-γ_a
    r_12 = S_a*τ_a*c_ab/N_b
    r_21 = S_b*τ_b*c_ab/N_a
    r_22 = S_b*τ_b*c_bb/N_b-γ_b

    tr = r_11+r_22
    det = r_11*r_22-r_12*r_21

    return (tr+sqrt(tr^2-4*det))/2
end
function calculate_inv_fit(mutant_strain::Real,resident_S,syspar::SystemParameters)  
    τ_a = τ_fun(mutant_strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_b = τ_fun(mutant_strain,syspar.z_b,syspar.σ²,syspar.τ_max)
    

    c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b, = syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.γ_a,syspar.γ_b,syspar.N_a,syspar.N_b

    S_a,S_b = resident_S[1],resident_S[2] #Susceptibles

    λ_max = calculate_inv_fit(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return λ_max
end

#Calculate Selection gradients (first order derivative of invasion fitness)
function calculate_select_grad(S_a,S_b,τ_a,τ_b,τ_prime_a,τ_prime_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    r_11 = S_a*τ_a*c_aa/N_a-γ_a
    r_12 = S_a*τ_a*c_ab/N_b
    r_21 = S_b*τ_b*c_ab/N_a
    r_22 = S_b*τ_b*c_bb/N_b-γ_b

    dr_11 = S_a*τ_prime_a*c_aa/N_a
    dr_12 = S_a*τ_prime_a*c_ab/N_b
    dr_21 = S_b*τ_prime_b*c_ab/N_a
    dr_22 = S_b*τ_prime_b*c_bb/N_b

    tr_R = r_11+r_22
    det_R =  r_11*r_22-r_12*r_21

    dtr_R = dr_11+dr_22
    ddet_R = dr_11*r_22+r_11*dr_22-dr_12*r_21-r_12*dr_21

    dλ_max = (dtr_R+(tr_R*dtr_R-2*ddet_R)/sqrt(tr_R^2-4*det_R))/2
    return dλ_max
end
function calculate_select_grad(strain::Real,S,syspar::SystemParameters)  
    τ_a = τ_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_b = τ_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)
    τ_prime_a = τ_prime_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_prime_b = τ_prime_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)

    c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b, = syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.γ_a,syspar.γ_b,syspar.N_a,syspar.N_b

    S_a,S_b = S[1],S[2] #Susceptibles

    dλ_max = calculate_select_grad(S_a,S_b,τ_a,τ_b,τ_prime_a,τ_prime_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return dλ_max
end

#Calculate second order derivative of invasion fitness
function calculate_sec_inv_fit(S_a,S_b,τ_a,τ_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    r_11 = S_a*τ_a*c_aa/N_a-γ_a
    r_12 = S_a*τ_a*c_ab/N_b
    r_21 = S_b*τ_b*c_ab/N_a
    r_22 = S_b*τ_b*c_bb/N_b-γ_b

    dr_11 = S_a*τ_prime_a*c_aa/N_a
    dr_12 = S_a*τ_prime_a*c_ab/N_b
    dr_21 = S_b*τ_prime_b*c_ab/N_a
    dr_22 = S_b*τ_prime_b*c_bb/N_b

    ddr_11 = S_a*τ_d_prime_a*c_aa/N_a
    ddr_12 = S_a*τ_d_prime_a*c_ab/N_b
    ddr_21 = S_b*τ_d_prime_b*c_ab/N_a
    ddr_22 = S_b*τ_d_prime_b*c_bb/N_b

    tr_R = r_11+r_22
    det_R =  r_11*r_22-r_12*r_21

    dtr_R = dr_11+dr_22
    ddet_R = dr_11*r_22+r_11*dr_22-dr_12*r_21-r_12*dr_21

    ddtr_R = ddr_11+ddr_22
    dddet_R = ddr_11*r_22+r_11*ddr_22+2*dr_11*dr_22-ddr_12*r_21-r_12*ddr_21-2*dr_12*dr_21

    ddλ_max = (ddtr_R+(tr_R*ddtr_R+dtr_R^2-2*dddet_R)/sqrt(tr_R^2-4*det_R)-(tr_R*dtr_R-2*ddet_R)^2/(tr_R^2-4*det_R)^(3/2))/2
    return ddλ_max
end
function calculate_sec_inv_fit(strain,S,syspar::SystemParameters)    
    τ_a = τ_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_b = τ_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)
    τ_prime_a = τ_prime_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_prime_b = τ_prime_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)
    τ_d_prime_a = τ_d_prime_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_d_prime_b = τ_d_prime_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)

    c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b, = syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.γ_a,syspar.γ_b,syspar.N_a,syspar.N_b
    
    S_a,S_b = S[1],S[2]

    ddλ_max = calculate_sec_inv_fit(S_a,S_b,τ_a,τ_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return ddλ_max
end

#Matrix (R) describing the dynamics of the two-species-one-pathogen system 
#dI = RI
function system_matrix(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)

    R = [S_a*τ_a*c_aa/N_a-γ_a S_a*τ_a*c_ab/N_b;S_b*τ_b*c_ab/N_a S_b*τ_b*c_bb/N_b-γ_b]
    return R
end
function system_matrix(strain,S,syspar::SystemParameters)
    τ_a = τ_fun(strain,syspar.z_a,syspar.σ²,syspar.τ_max)
    τ_b = τ_fun(strain,syspar.z_b,syspar.σ²,syspar.τ_max)

    c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b, = syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.γ_a,syspar.γ_b,syspar.N_a,syspar.N_b

    S_a,S_b = S[1],S[2]

    R = system_matrix(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return R
end

#Check satablility of system matrix of the two-species-one-pathogen system
#True if system R0>1 for strain
function check_system_matrix(strain,S,syspar::SystemParameters)    
    R = system_matrix(strain,S,syspar)
    eigs = eigvals(R)
    return any(real.(eigs).>0)
end
function check_system_matrix(strain,syspar::SystemParameters)        
    return check_system_matrix(strain,(syspar.N_a,syspar.N_b),syspar)    
end

#System dynamics for the two-species-one-pathogen system
function epievodyn_simple_one_strain!(du,u,p,t)
    strain,syspar = p    

    S = u[1:2,1]
    I = u[1:2,2]

    R = system_matrix(strain,S,syspar) 
        
    du[1:2,1] = -R*I
    du[1:2,2] = R*I
end

#Function to find the equilibrium point for the two-species-one-pathogen system (OLD To be replaced)
function find_eq(u₀,tspan,p)
    t_end = tspan[2] 
    prob = ODEProblem(epievodyn_simple_one_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    S_a_eq,S_b_eq = sol(t_end)[1,1],sol(t_end)[2,1]
    I_a_eq,I_b_eq = sol(t_end)[1,2],sol(t_end)[2,2]

    eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)

    return eq_point
end

#-Start: Code and functions to find equilibrium points for the two-species-one-strain system

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
#The roots of f(x) are equilibrium point candidates
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

#System dynamics for the two-species-two-pathogen system
function epievodyn_simple_2_strain!(du,u,p,t)
    strain_1,strain_2,syspar = p
  
    S = u[1:2,1]
    I_1 = u[1:2,2]
    I_2 = u[1:2,3]

    R_1 = system_matrix(strain_1,S,syspar) 
    R_2 = system_matrix(strain_2,S,syspar) 

    du[1:2,1] = -R_1*I_1-R_2*I_2
    du[1:2,2] = R_1*I_1
    du[1:2,3] = R_2*I_2 
end

#Function to find the equilibrium point for the two-species-two-pathogen system (NEED TO UPDATE)
function find_eq_2_strain(u₀,tspan,p)
    t_end = tspan[2] 
    prob = ODEProblem(epievodyn_simple_2_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    S_a_eq,S_b_eq = sol(t_end)[1,1],sol(t_end)[2,1]
    I_a_eq,I_b_eq = sol(t_end)[1,2]+sol(t_end)[1,3],sol(t_end)[2,2]+sol(t_end)[2,3]

    eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)

    return eq_point
end

function create_PIP(syspar::SystemParameters,option)
    #z_start,z_end,Nₚ,u₀,tspan = option 
    z_start,z_end,Nₚ = option   

    strategies = range(z_start,stop=z_end,length=Nₚ) 
    
    r_m_mat = zeros(Nₚ,Nₚ) #PIP plot
    
    for col_idx in 1:1:Nₚ #Resident strain
        resident_strain = strategies[col_idx]
               
        if check_system_matrix(resident_strain,(syspar.N_a,syspar.N_b),syspar) 
            #eq_point = find_eq(u₀,tspan,(resident_strain,syspar))
            eq_point = find_eq(resident_strain,syspar)
            S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
            for i in 1:1:Nₚ #Mutant strain strain
                mutant_strain = strategies[i]  
                                
                if check_system_matrix(mutant_strain,(S_a_eq,S_b_eq),syspar)
                    r_m_mat[i,col_idx]=1
                end
            end
        else
            r_m_mat[:,col_idx].=0.5            
        end        
    end
    return r_m_mat
end

function draw_PIP(fig_name,strategies,r_m_mat,z_start,z_end,z_a,z_b)
    heatmap(strategies,strategies,r_m_mat)
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([z_a,z_a],[z_start,z_end],c=:red)
    plot!([z_b,z_b],[z_start,z_end],xlims=[z_start,z_end],legend=false,c=:blue)
    savefig(fig_name) 
    return nothing
end
function draw_PIP(fig_name,syspar::SystemParameters,option)
    #z_start,z_end,Nₚ,u₀,tspan = option  
    z_start,z_end,Nₚ = option   

    strategies = range(z_start,stop=z_end,length=Nₚ)
    
    r_m_mat = create_PIP(syspar,option)

    draw_PIP(fig_name*".svg",strategies,r_m_mat,z_start,z_end,syspar.z_a,syspar.z_b)
    
    return nothing
end

#Create map of of coexisting strains from PIP
function create_coex_region(pip::Matrix{T}) where{T}
    n_row,n_col = size(pip)
    n_row==n_col||error("n_row!=n_col")
    coex_region = zeros(T,n_row,n_col)
    for row in 1:n_row
        for col in 1:n_col 
            r_1 = pip[row,col]
            r_2 = pip[col,row]
            if r_1≈0.5
                coex_region[row,col] = 0.5
            else
                if r_1≈1.0&&r_2≈1.0
                    coex_region[row,col]=1
                else
                    coex_region[row,col]=0
                end
            end
        end
    end
    return coex_region
end
function create_coex_region(syspar::SystemParameters,option)
    r_m_mat = create_PIP(syspar,option)
    coex_region = create_coex_region(r_m_mat)
    return coex_region
end

#Draw region of coexisting strains
function draw_coex_region(fig_name,syspar::SystemParameters,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    μ_a,μ_b = syspar.z_a,syspar.z_b

    coex_region = create_coex_region(syspar,option)

    z_vec = range(z_start,stop=z_end,length=Nₚ) 
    heatmap(z_vec,z_vec,coex_region)   
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    plot!([z_start,z_end],[μ_a,μ_a],c=:red)
    plot!([z_start,z_end],[μ_b,μ_b],c=:blue)
    savefig(fig_name*".svg")

    return nothing
end

function create_TEP(syspar::SystemParameters,option,coex_region::Matrix{T};ϵ::Real=1e-6) where{T}
    z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain = option
        
    strategies = range(z_start,stop=z_end,length=Nₚ)

    (n_row,n_col) = size(coex_region)
    n_row==n_col==Nₚ||error("n_row!=n_col!=Nₚ")    

    tep = zeros(T,n_row,n_col)

    for row in 1:n_row
        for col in 1:n_col            
            if coex_region[row,col]≈0.5
                tep[row,col] = 0.5
            elseif coex_region[row,col]≈0
                tep[row,col] = 0
            else
                r_1 = strategies[row] #Resident strain 1
                r_2 = strategies[col] #Resident strain 2                
                p_system = r_1,r_2,syspar
                eq = find_eq_2_strain(u₀_2_strain,tspan,p_system)
            
                dλ_max_1 = calculate_select_grad(r_1,(eq.S_a,eq.S_b),syspar)
                dλ_max_2 = calculate_select_grad(r_2,(eq.S_a,eq.S_b),syspar)

                ddλ_max_1 = calculate_sec_inv_fit(r_1,(eq.S_a,eq.S_b),syspar)
                ddλ_max_2 = calculate_sec_inv_fit(r_2,(eq.S_a,eq.S_b),syspar)

                if abs(dλ_max_1)<ϵ 
                    if  ddλ_max_1<0 #fitness maximum
                        tep[row,col]=0.7
                    else
                        tep[row,col]=0.3
                    end
                elseif abs(dλ_max_2)<ϵ
                    if  ddλ_max_2<0 #fitness maximum
                        tep[row,col]=0.7
                    else
                        tep[row,col]=0.3
                    end                   
                else
                    tep[row,col]=1
                end
            end
        end
    end

    return tep
end

function calc_invasion_cone(r_1,r_2,syspar::SystemParameters,option)
    z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain = option
    
    p_system = r_1,r_2,syspar   
    eq = find_eq_2_strain(u₀_2_strain,tspan,p_system)
       
    dλ_max_r1 = calculate_select_grad(r_1,(eq.S_a,eq.S_b),syspar)  
    dλ_max_r2 = calculate_select_grad(r_2,(eq.S_a,eq.S_b),syspar) 
    
    return (dλ_max_r1,dλ_max_r2)
end

#Find equilibrium point and calculate selection gradiant at this point
function find_selectgrad(strain::Real,syspar::SystemParameters,option)
    u₀,tspan = option
    p = strain,syspar 
    
    eq_point = find_eq(u₀,tspan,p)
    
    S = (eq_point.S_a,eq_point.S_b)
    dλ_max = calculate_select_grad(strain,S,syspar)

    return dλ_max
end

#Find equilibrium point and calculate second order derivative of invasion fitness at this point
function find_sec_inv_fit(strain::Real,syspar::SystemParameters,option)
    u₀,tspan = option
    p = strain,syspar 
    
    eq_point = find_eq(u₀,tspan,p)
    
    S = (eq_point.S_a,eq_point.S_b)
    
    ddλ_max = calculate_sec_inv_fit(strain, S, syspar)

    return ddλ_max
end

#Check if singular strategy strain is convergence stable
function check_convergence(strain::Real,syspar::SystemParameters,option;h::Real=10^-3)
    
    dλ_max_1 = find_selectgrad(strain+h,syspar,option)
    dλ_max_2 = find_selectgrad(strain-h,syspar,option)

    return (dλ_max_1-dλ_max_2)/(2*h)<0
end

#Find bracketing intervals (sign change intervals) in the range  defined [xstart,xend]
function find_bracketing_interval(fun::Function,xstart::Real,xend::Real,length::Integer)
    xend>xstart||error("xend≤xstart")
    length>=2||error("length≤2")

    Δx = (xend-xstart)/(length-1)             
    
    bracketing = Tuple{Real,Real}[]
    left_val = xstart
    left_sign = fun(left_val)<0  

    for i = 1:length-1 
        right_val = xstart+Δx*i
        right_sign = fun(right_val)<0
        if right_sign!=left_sign
            
            push!(bracketing,(left_val,right_val)) 

            left_sign = right_sign
            left_val = right_val
        end
    end

    if isempty(bracketing)
        println("No sign change detected in the given intervall")
        return nothing  
    end

    return bracketing    
end

function find_singular_strategies(syspar::SystemParameters,option)
    z_start,z_end,Nₚ,u₀,tspan = option    
    option_temp = u₀,tspan
    Nₚ>2||error("Nₚ≤2")      
    b_vals = find_bracketing_interval(x->find_selectgrad(x,syspar,option_temp),z_start,z_end,Nₚ)
    isnothing(b_vals)==false||error("No sign change detected in the given intervall")
    z_0s = zeros(length(b_vals))    
    for (i,val) in enumerate(b_vals)        
        z_0s[i] = bisection(x->find_selectgrad(x,syspar,option_temp), val[1],val[2])             
    end
    return z_0s
end

function singular_strategies(syspar::SystemParameters,option)
    z_start,z_end,Nₚ,u₀,tspan = option 

    option_in = u₀,tspan

    Δz = (z_end-z_start)/(Nₚ-1)
      
    z0s = find_singular_strategies(syspar,option)  

    sing_strats = SingularStrat[]

    for z0 in z0s         
        conv_check = check_convergence(z0,syspar,option_in;h=Δz)    
        evo_check = find_sec_inv_fit(z0,syspar,option_in)
        push!(sing_strats,SingularStrat(z0,conv_check,evo_check<0))
    end

    return sing_strats
end

function draw_PIP_n_coex_region(fig_name::AbstractString,syspar::SystemParameters)  
    z_start = syspar.z_a-0.02
    z_start>0||error("z_start<0")
    z_end = syspar.z_b+0.02
    
    #Initial values
    Nₚ = 1000
    N_a = syspar.N_a    
    N_a>2||error("N_a<3")
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
     
    N_b = syspar.N_b   
    N_b>2||error("N_b<3")
    I_b₀ = 1
    S₀_b = N_b-I_b₀  

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]  

    option = z_start,z_end,Nₚ,u₀,tspan 

    draw_PIP(fig_name*"PIP",syspar,option)    
    draw_coex_region(fig_name*"Coex_reg",syspar,option)  

    return nothing
end

function check_R₀_region(syspar::SystemParameters,z_start::Real,z_end::Real;N::Integer=1000)    
    S₀ = [syspar.N_a,syspar.N_b]
    strains=range(z_start,stop=z_end,length=N)

    check_system_matrix(z_start,S₀,syspar)&&check_system_matrix(z_end,S₀,syspar)||error("R₀ for one of the species is < 1.0")

    matrix_instable = [check_system_matrix(strain,S₀,syspar) for strain in strains]
    
    intervals = Tuple{Real,Real,Bool}[]

    start_interval_idx = 1
    interval_type = matrix_instable[1]

    for i in eachindex(strains)
        if interval_type!=matrix_instable[i]
            start_interval_idx!=i-1||error("interval error: try to increas N")
            push!(intervals,(strains[start_interval_idx],strains[i-1],interval_type))
            start_interval_idx  = i
            interval_type = matrix_instable[i]
        end
    end

    push!(intervals,(strains[start_interval_idx],strains[end],interval_type))

    return intervals
end

function get_singular_strategies(syspar::SystemParameters,
    z_start::Real,
    z_end::Real;
    Nₚ::Integer = 3000,
    t_end::Real = 8000)
    
    R₀_intervals = check_R₀_region(syspar,z_start,z_end;N=Nₚ)  
     
    N_a = syspar.N_a    
    N_a>2||error("N_a<3")
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
    
    N_b = syspar.N_b   
    N_b>2||error("N_b<3")
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    strat = SingularStrat[]

    for interval in R₀_intervals 
        if last(interval)
            z_start_temp,z_end_temp = interval[1],interval[2]
            option = z_start_temp,z_end_temp,Nₚ,u₀,tspan 
        
            strat_temp = singular_strategies(syspar,option)
            append!(strat,strat_temp)
        end
    end
    return (strat,any(last.(R₀_intervals).==false))
end