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

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
τ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
τ_prime_fun(x,μ,σ²,amplitude) = -amplitude*exp(-(x-μ)^2/(σ²*2))*2*(x-μ)/(σ²*2)
τ_d_prime_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))*4*(x-μ)^2/(σ²*2)^2-amplitude*exp(-(x-μ)^2/(σ²*2))*2/(σ²*2)

#Calculate Selection gradients
function calculate_select_grad(S_a,S_b,τ_a,τ_b,τ_prime_a,τ_prime_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    r_11 = S_a*τ_a*c_aa/N_a-γ_a
    r_12 = S_a*τ_a*c_ab/N_b
    r_21 = S_b*τ_b*c_ab/N_a
    r_22 = S_b*τ_b*c_bb/N_b-γ_b

    dr_11 = S_a*τ_prime_a*c_aa/N_a
    dr_12 = S_a*τ_prime_a*c_ab/N_b
    dr_21 = S_b*τ_prime_b*c_ab/N_a
    dr_22 = S_b*τ_prime_b*c_bb/N_b

    R = [r_11 r_12;r_21 r_22]
    dR = [dr_11 dr_12;dr_21 dr_22]

    R_e_val,R_e_vec_R = eigen(R)
    R_e_val,R_e_vec_L = eigen(transpose(R))

    R_e_vec_L_T = transpose(R_e_vec_L[:,2])

    dλ_max_alt = R_e_vec_L_T*dR*R_e_vec_R[:,2]/(R_e_vec_L_T*R_e_vec_R[:,2])

    tr_R = r_11+r_22
    det_R =  r_11*r_22-r_12*r_21

    dtr_R = dr_11+dr_12
    ddet_R = dr_11*r_22+r_11*dr_22-dr_12*r_21-r_12*dr_21

    dλ_max = (dtr_R+(tr_R*dtr_R-2*ddet_R)/sqrt(tr_R^2-4*det_R))/2
    return dλ_max
end
function calculate_select_grad(p,S)    
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b = p
    S_a,S_b = S[1],S[2]
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

    dtr_R = dr_11+dr_12
    ddet_R = dr_11*r_22+r_11*dr_22-dr_12*r_21-r_12*dr_21

    ddtr_R = ddr_11+ddr_12
    dddet_R = ddr_11*r_22+r_11*ddr_22+2*dr_11*dr_22-ddr_12*r_21-r_12*ddr_21-2*dr_12*dr_21

    ddλ_max = (ddtr_R+(tr_R*ddtr_R+dtr_R^2-2*dddet_R)/sqrt(tr_R^2-4*det_R)-(tr_R*dtr_R-2*ddet_R)^2/(tr_R^2-4*det_R)^(3/2))/2
    return ddλ_max
end
function calculate_sec_inv_fit(p,S)    
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b = p
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
function system_matrix(p,S)
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b = p
    S_a,S_b = S[1],S[2]
    R = system_matrix(S_a,S_b,τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b)
    return R
end

#Check satablility of system matrix of the two-species-one-pathogen system 
function check_system_matrix(p,S)    
    R = system_matrix(p,S)
    eigs = eigvals(R)
    return any(real.(eigs).>0)
end

#System dynamics for the two-species-one-pathogen system
function epievodyn_simple_one_strain!(du,u,p,t)
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b = p

    S = u[1:2,1]
    I = u[1:2,2]

    R = system_matrix(S[1],S[2],τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b) 
        
    du[1:2,1] = -R*I
    du[1:2,2] = R*I
end

#Function to find the equilibrium point for the two-species-one-pathogen system
function find_eq(u₀,tspan,p)
    t_end = tspan[2] 
    prob = ODEProblem(epievodyn_simple_one_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    S_a_eq,S_b_eq = sol(t_end)[1,1],sol(t_end)[2,1]
    I_a_eq,I_b_eq = sol(t_end)[1,2],sol(t_end)[2,2]

    eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)

    return eq_point
end

#System dynamics for the two-species-two-pathogen system
function epievodyn_simple_two_strain!(du,u,p,t)
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,strain_1,strain_2 = p
  
    S = u[1:2,1]
    I_1 = u[1:2,2]
    I_2 = u[1:2,3]

    R_1 = system_matrix(S[1],S[2],τ_a(strain_1),τ_b(strain_1),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b) 
    R_2 = system_matrix(S[1],S[2],τ_a(strain_2),τ_b(strain_2),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b) 

    du[1:2,1] = -R_1*I_1-R_2*I_2
    du[1:2,2] = R_1*I_1
    du[1:2,3] = R_2*I_2 
end

#Function to find the equilibrium point for the two-species-two-pathogen system
function find_eq_2_strain(u₀,tspan,p)
    t_end = tspan[2] 
    prob = ODEProblem(epievodyn_simple_two_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    S_a_eq,S_b_eq = sol(t_end)[1,1],sol(t_end)[2,1]
    I_a_eq,I_b_eq = sol(t_end)[1,2]+sol(t_end)[1,3],sol(t_end)[2,2]+sol(t_end)[2,3]

    eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)

    return eq_point
end


#Calcualte the R0 (reproduction number disease free state) for the two-species-one-pathogen system 
function calc_R0(p)
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b = p
   
    K_L = [τ_a*c_aa/γ_a τ_a*c_ab*N_a/(N_b*γ_b);τ_b*c_ab*N_b/(N_a*γ_a) τ_b*c_bb/γ_b]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end

#Calcualte the R0 for a rare mutant pathogen in a environment set by a recident pathogen 
function calc_R0_m(p,S_a_eq,S_b_eq)
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b = p
    
    K_L = [τ_a*c_aa*S_a_eq/(γ_a*N_a) τ_a*c_ab*S_a_eq/(N_b*γ_b);τ_b*c_ab*S_b_eq/(N_a*γ_a) τ_b*c_bb*S_b_eq/(γ_b*N_b)]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end
 
function create_PIP(p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,z_a,z_b = p_in

    strategies = range(z_start,stop=z_end,length=Nₚ) 
    
    r_m_mat = zeros(Nₚ,Nₚ) #PIP plot
    r_m_mat_alt = zeros(Nₚ,Nₚ) #PIP plot

    for col_idx in 1:1:Nₚ #Resident strain
        resident_strain = strategies[col_idx]

        p_resident = τ_a(resident_strain),τ_b(resident_strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b
        R0_resident = calc_R0(p_resident)
        if R0_resident>1.0
            eq_point = find_eq(u₀,tspan,p_resident)
            S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
            for i in 1:1:Nₚ #Mutant strain strain
                mutant_strain = strategies[i]  
                p_mutant = τ_a(mutant_strain),τ_b(mutant_strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b

                if calc_R0_m(p_mutant,S_a_eq,S_b_eq)>1.0
                    r_m_mat[i,col_idx]=1
                end
                if check_system_matrix(p_mutant,(S_a_eq,S_b_eq))
                    r_m_mat_alt[i,col_idx]=1
                end
            end
        else
            r_m_mat[:,col_idx].=0.5
            r_m_mat_alt[:,col_idx].=0.5
        end        
    end
    return r_m_mat,r_m_mat_alt
end

function draw_PIP(fig_name,strategies,r_m_mat,z_start,z_end,z_a,z_b)
    heatmap(strategies,strategies,r_m_mat)
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([z_a,z_a],[z_start,z_end],c=:red)
    plot!([z_b,z_b],[z_start,z_end],legend=false,c=:blue)
    savefig(fig_name) 
    return nothing
end
function draw_PIP(fig_name,p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,z_a,z_b = p_in

    strategies = range(z_start,stop=z_end,length=Nₚ)
    
    r_m_mat,r_m_mat_alt = create_PIP(p_in,option)

    draw_PIP(fig_name*".svg",strategies,r_m_mat,z_start,z_end,z_a,z_b)
    draw_PIP(fig_name*"_alt.svg",strategies,r_m_mat_alt,z_start,z_end,z_a,z_b)
    
    return nothing
end

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

function create_coex_region(p_in,option)
    r_m_mat,r_m_mat_alt = create_PIP(p_in,option)
    coex_region = create_coex_region(r_m_mat)
    return coex_region
end

function create_TEP(p_in,option,coex_region::Matrix{T};ϵ::Real=1e-6) where{T}
    z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,z_a,z_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b = p_in
    
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
                p_system = τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,r_1,r_2
                eq = find_eq_2_strain(u₀_2_strain,tspan,p_system)

                p_grad_1 = τ_a(r_1),τ_b(r_1),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_1),τ_prime_b(r_1)
                p_grad_2 = τ_a(r_2),τ_b(r_2),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_2),τ_prime_b(r_2)

                p_sec_deriv_1 = τ_a(r_1),τ_b(r_1),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_1),τ_prime_b(r_1),τ_d_prime_a(r_1),τ_d_prime_b(r_1)
                p_sec_deriv_2 = τ_a(r_2),τ_b(r_2),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_2),τ_prime_b(r_2),τ_d_prime_a(r_2),τ_d_prime_b(r_2)

                dλ_max_1 = calculate_select_grad(p_grad_1,(eq.S_a,eq.S_b))
                dλ_max_2 = calculate_select_grad(p_grad_2,(eq.S_a,eq.S_b))

                ddλ_max_1 = calculate_sec_inv_fit(p_sec_deriv_1,(eq.S_a,eq.S_b))
                ddλ_max_2 = calculate_sec_inv_fit(p_sec_deriv_2,(eq.S_a,eq.S_b))

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

function calc_invasion_cone(r_1,r_2,p_in,option)
    z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,z_a,z_b,τ_prime_a,τ_prime_b = p_in

    p_system = τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,r_1,r_2    
    eq = find_eq_2_strain(u₀_2_strain,tspan,p_system)
   
    p_grad_r1 = τ_a(r_1),τ_b(r_1),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_1),τ_prime_b(r_1)
    dλ_max_r1 =  calculate_select_grad(p_grad_r1,(eq.S_a,eq.S_b))    
    
    p_grad_r2 = τ_a(r_2),τ_b(r_2),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(r_2),τ_prime_b(r_2)
    dλ_max_r2 =  calculate_select_grad(p_grad_r2,(eq.S_a,eq.S_b)) 
    
    return (dλ_max_r1,dλ_max_r2)
end