#Simple SIS-model example with two species and one strategy (z:=[0,1])

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
γ_fun_prim(x,μ,σ²,amplitude) = -(x-μ)/(σ²*2)*amplitude*exp(-(x-μ)^2/(σ²*2))

function system_matrix(p,S)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    S_a,S_b = S[1],S[2]
    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)
    r = [S_a*γ_aa*K_aa/N_a-c S_a*γ_aa*K_ab/N_b;S_b*γ_bb*K_ab/N_a S_b*γ_bb*K_bb/N_b-c]
    return r
end
function system_matrix(p,mutant_strain,S_eq)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    S_a,S_b = S_eq[1],S_eq[2]
    γ_aa = γ_a(mutant_strain)
    γ_bb = γ_b(mutant_strain)
    r = [S_a*γ_aa*K_aa/N_a-c S_a*γ_aa*K_ab/N_b;S_b*γ_bb*K_ab/N_a S_b*γ_bb*K_bb/N_b-c]
    return r
end

function check_system_matrix(p,S)
    r = system_matrix(p,S)
    eigs = eigvals(r)
    return any(real.(eigs).>0)
end
function check_system_matrix(p,mutant_strain,S_eq)
    r = system_matrix(p,mutant_strain,S_eq)
    eigs = eigvals(r)
    return any(real.(eigs).>0)
end

function epievodyn_simple_one_strain!(du,u,p,t)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p

    S_a,S_b = u[1,1],u[1,2]
    I_a, I_b = u[2,1],u[2,2]
       
    temp_a = S_a*γ_a(strategy)*(I_a*K_aa/N_a+I_b*K_ab/N_b)
    temp_b = S_b*γ_b(strategy)*(I_a*K_ab/N_a+I_b*K_bb/N_b)
        
    du[1,1] = c*I_a-temp_a 
    du[1,2] = c*I_b-temp_b

    du[2,1] = -du[1,1]
    du[2,2] = -du[1,2]

    return du
end

mutable struct EqPoint{T<:Real}
    S_a::T
    S_b::T
    I_a::T
    I_b::T  
end

function find_eq(u₀,tspan,p)
    t_end = tspan[2] 
    prob = ODEProblem(epievodyn_simple_one_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE
    
    du = [0.0 0.0;0.0 0.0]
    @show u = sol(t_end)
    @show epievodyn_simple_one_strain!(du,sol(t_end),p,t_end)

    @show S_a_eq,S_b_eq = sol(t_end)[1,1],sol(t_end)[1,2]
    @show I_a_eq,I_b_eq = sol(t_end)[2,1],sol(t_end)[2,2]

    eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)

    return eq_point
end

function system_res!(F,I,p)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p

    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)

    I_a = I[1]
    I_b = I[2]

    S_a = N_a-I_a
    S_b = N_b-I_b

    temp_a = S_a*γ_aa*(I_a*K_aa/N_a+I_b*K_ab/N_b)
    temp_b = S_b*γ_bb*(I_a*K_ab/N_a+I_b*K_bb/N_b)      
       
    F[1] = temp_a-c*I_a
    F[2] = temp_b-c*I_b
    
end
function system_res_J!(J,I,p)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p

    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)

    I_a = I[1]
    I_b = I[2]

    S_a = N_a-I_a
    S_b = N_b-I_b  
    
    #dI_a=-I_a^2*γ_aa*K_aa/N_a-I_a*I_b*γ_aa*K_ab/N_b+I_a*(γ_aa*K_aa-c)+I_b*N_a*γ_aa*K_ab/N_b
    #dI_a=-I_b^2*γ_bb*K_bb/N_b-I_a*I_b*γ_bb*K_ab/N_a+I_b*(γ_bb*K_bb-c)+I_a*N_b*γ_bb*K_ab/N_a

    J[1,1] = -2*I_a*γ_aa*K_aa/N_a-I_b*γ_aa*K_ab/N_b+(γ_aa*K_aa-c)
    J[1,2] = -I_a*γ_aa*K_ab/N_b+N_a*γ_aa*K_ab/N_b
    J[2,1] = -I_b*γ_bb*K_ab/N_a+N_b*γ_bb*K_ab/N_a
    J[2,2] = -2*I_b*γ_bb*K_bb/N_b-I_a*γ_bb*K_ab/N_a+(γ_bb*K_bb-c)

end

#Find equilibrium point
function find_eq(I₀,p)     

   f!(F,I) = system_res!(F,I,p)  
   j!(F,I) = system_res_J!(F,I,p) 
   
    ans = nlsolve(f!, j!, I₀)

    #eq_point = EqPoint(S_a_eq,S_b_eq,I_a_eq,I_b_eq)
    return ans
end

function find_new_strategy(u₀,tspan,p,γ_prim)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_a_prim,γ_b_prim = γ_prim

    eq_point = find_eq(u₀,tspan,p)
        
    S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
    I_a_eq,I_b_eq = eq_point.I_a,eq_point.I_b

    I = I_a_eq+I_b_eq
    q_a = I_a_eq/I
    q_b = I_b_eq/I

    K_N_mean_a =  K_aa/N_a*q_a+K_ab/N_b*q_b
    K_N_mean_b = K_ab/N_a*q_a+K_bb/N_b*q_b     
    
    r_m_prim_fun(z) = γ_a_prim(z)*S_a_eq*K_N_mean_a+γ_b_prim(z)*S_b_eq*K_N_mean_b

    r_m_prim = r_m_prim_fun(strategy)
    u₀_new = [S_a_eq S_b_eq;I_a_eq I_b_eq]

    return (strategy*(1+sign(r_m_prim)*0.001),u₀_new)
end

function calc_R0(p)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)
    K_L = [γ_aa*K_aa/c γ_aa*K_ab*N_a/(N_b*c);γ_bb*K_ab*N_b/(N_a*c) γ_bb*K_bb/c]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end

function draw_vector_field(p,I_a_start,I_a_end,I_a_N,I_b_start,I_b_end,I_b_N)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_aa = γ_a(strategy)
    γ_bb = γ_b(strategy)

    I_a_vec = ones(I_a_N*I_b_N)
    I_b_vec = ones(I_a_N*I_b_N)
    dI_a_vec = ones(I_a_N*I_b_N)
    dI_b_vec = ones(I_a_N*I_b_N)

    for (i,I_a) in enumerate(range(I_a_start,stop=I_a_end,length=I_a_N))
        for (j,I_b) in enumerate(range(I_b_start,stop=I_b_end,length=I_b_N))
            idx = (i-1)*I_b_N+j

            I_a_vec[idx] = I_a
            I_b_vec[idx] = I_b

            S_a = N_a-I_a
            S_b = N_b-I_b

            temp_a = S_a*γ_aa*(I_a*K_aa/N_a+I_b*K_ab/N_b)
            temp_b = S_b*γ_bb*(I_a*K_ab/N_a+I_b*K_bb/N_b)
                
            dI_a = temp_a-c*I_a
            dI_b = temp_b-c*I_b

            dI_a_vec[idx] = dI_a
            dI_b_vec[idx] = dI_b
        end
    end

    quiver(I_a_vec,I_b_vec,quiver=(dI_a_vec,dI_b_vec))
end

function calc_R0_m(p,S_a_eq,S_b_eq,strategy_m)
    c,γ_a,γ_b,strategy_r,K_aa,K_bb,K_ab,N_a,N_b = p
    γ_aa = γ_a(strategy_m)
    γ_bb = γ_b(strategy_m)
    K_L = [γ_aa*K_aa*S_a_eq/(c*N_a) γ_aa*K_ab*S_a_eq/(N_b*c);γ_bb*K_ab*S_b_eq/(N_a*c) γ_bb*K_bb*S_b_eq/(c*N_b)]

    return abs((tr(K_L)+sqrt(tr(K_L)^2-4*det(K_L)))/2)
end

function draw_curvature(p,u₀,tspan,delta;n_point::Integer=100)
    c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b = p

    eq_point = find_eq(u₀,tspan,p)

    S_a_eq,S_b_eq = eq_point.S_a,eq_point.S_b
    I_a_eq,I_b_eq = eq_point.I_a,eq_point.I_b

    I = I_a_eq+I_b_eq
    q_a = I_a_eq/I
    q_b = I_b_eq/I

    K_N_mean_a =  K_aa/N_a*q_a+K_ab/N_b*q_b
    K_N_mean_b = K_ab/N_a*q_a+K_bb/N_b*q_b 

    r_m(z) = γ_a(z)*S_a_eq*K_N_mean_a+γ_b(z)*S_b_eq*K_N_mean_b-c 

    z_vec = range(strategy-delta,stop=strategy+delta,length=n_point)

    return (r_m.(z_vec),z_vec)
end

function run()
    Nₚ = 200
    z_vec = range(0.0,stop=1.0,length=Nₚ+1) 

    #Parameters
    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3 

    c = 0.1    

    #Initial values    
    N_a = 10^3  
    I_a₀ = 2  
    S₀_a = N_a-I_a₀ 
    N_b = 10^3     
    I_b₀ = 0
    S₀_b = N_b-I_b₀ 

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.6
    μ_b,σ²_b,amplitude_b = 0.35,0.0025,0.6

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b)
    γ_a_prim(z) = γ_fun_prim(z,μ_a,σ²_a,amplitude_a)
    γ_b_prim(z) = γ_fun_prim(z,μ_b,σ²_b,amplitude_b)

    γ_prim = γ_a_prim,γ_b_prim

    strategy = 0.2

    p = c,γ_a,γ_b,strategy,K_aa,K_bb,K_ab,N_a,N_b

    t_start = 0
    t_end = 900
    tspan = (t_start,t_end)

    u₀ = [S₀_a S₀_b;I_a₀ I_b₀]

    prob = ODEProblem(epievodyn_simple_one_strain!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    @show strategy_new,u₀_new = find_new_strategy(u₀,tspan,p,γ_prim)

    for i = 1:2551
        p_new = c,γ_a,γ_b,strategy_new,K_aa,K_bb,K_ab,N_a,N_b
        @show strategy_new,u₀_new = find_new_strategy(u₀_new,tspan,p_new,γ_prim)
    end     
    
    
    @show S_a_eq,S_b_eq = u₀_new[1,1],u₀_new[1,2]
    @show I_a_eq,I_b_eq = u₀_new[2,1],u₀_new[2,2]

    @show I = I_a_eq+I_b_eq
    @show q_a = I_a_eq/I
    @show q_b = I_b_eq/I

    @show K_N_mean_a =  K_aa/N_a*q_a+K_ab/N_b*q_b
    @show K_N_mean_b = K_ab/N_a*q_a+K_bb/N_b*q_b 
    
    r_m(z) = γ_a(z)*S_a_eq*K_N_mean_a+γ_b(z)*S_b_eq*K_N_mean_b-c    

    strategy_m = range(strategy_new*0.99,stop=strategy_new*1.01,length=100)  
    plot(strategy_m,r_m.(strategy_m))
    plot!([strategy_new,strategy_new],[minimum(r_m.(strategy_m)),maximum(r_m.(strategy_m))])   

    
end

#run()