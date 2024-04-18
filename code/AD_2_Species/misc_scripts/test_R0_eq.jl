function R₀_fun(z::Real,syspar::SystemParameters)

    x_a,x_b,c_aa,c_bb,c_ab,N_a,N_b,γ_a,γ_b,σ²,τ_max = syspar.z_a,syspar.z_b,syspar.c_aa,syspar.c_bb,syspar.c_ab,syspar.N_a,syspar.N_b,syspar.γ_a,syspar.γ_b,syspar.σ²,syspar.τ_max

    R₀_a_max = τ_max*c_aa/γ_a
    R₀_b_max = τ_max*c_bb/γ_b
    R₀_ab_max = τ_max*c_ab/γ_b
    R₀_ba_max = τ_max*c_ab/γ_a

    R₀_a = R₀_a_max*τ_fun(x_a,z,σ²) 
    R₀_b = R₀_b_max*τ_fun(x_b,z,σ²) 
    R₀_ab = R₀_ab_max*τ_fun(x_a,z,σ²) 
    R₀_ba = R₀_ba_max*τ_fun(x_b,z,σ²) 

    R₀_mean =  (R₀_a+R₀_b)/2
    ΔR₀ =  (R₀_a-R₀_b)/2

    return R₀_mean+sqrt(ΔR₀^2+R₀_ab*R₀_ba) 
end

function test_R0()
    γ = 0.1 
    z_start = 0.18
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.33
    Nₚ = 1500    

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.2
        
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    t_end = 300
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]  

    syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)

    @show R₀_a_max = amplitude*c_aa/γ
    @show R₀_b_max = amplitude*c_bb/γ
    @show R₀_ab_max = amplitude*c_ab/γ

    Δz = range(0.0,stop=1.0,length=200)
    Δx = μ_b-μ_a
    z = [Δx*Δ+μ_a for Δ in Δz]    
    plot(Δz,R₀_fun.(z,Ref(syspar)))

    options = (0.2,syspar)
    prob = ODEProblem(epievodyn_simple_one_strain!,u₀,tspan,options) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE
    plot(sol)

    @show eq = find_eq(u₀,tspan,options)
    check_system_matrix(0.2,[eq.S_a,eq.S_b],syspar)  

    probSS = SteadyStateProblem(epievodyn_simple_one_strain!,u₀,options)
    solss = solve(probSS, DynamicSS(Rodas5(),abstol=1e-7,reltol=1e-5))  
    @show solss
    @show u = [solss[1] solss[3];solss[2] solss[4]]
    du = zeros(2,2)
    epievodyn_simple_one_strain!(du,u,options,0.0)

    du
end 

test_R0()