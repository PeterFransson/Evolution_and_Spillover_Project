#Simple SIS-model example with two species and one strategy (z:=[0,1])

#Simpson's 1/3 rule
simpson13(f,a,b) = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6

function cont2disc(f,disc_samples)
    n = length(disc_samples)
    areas = simpson13.(Ref(f),disc_samples[1:n-1],disc_samples[2:n])
    area_tot = sum(areas)
    return areas/area_tot
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))
m(z,z_vec,μₘ,σ²ₘ) = normal_d(z,μₘ,σ²ₘ)/sum(normal_d.(z_vec,Ref(μₘ),Ref(σ²ₘ)))

function ρ_fun(zᵢ,zⱼ,μ,σ²,amplitude)
    if γ_fun(zᵢ,μ,σ²,amplitude)<γ_fun(zⱼ,μ,σ²,amplitude)
        return 0;
    else
        return 1; 
    end
end

function epievodyn_simple!(du,u,p,t)
    c,γ_a,γ_b,strategies,Nₚ,K_aa,K_bb,K_ab,N_a,N_b = p

    S_a,S_b = u[1,1],u[1,2]
    I_a_vec = u[2:Nₚ+1,1]
    I_b_vec = u[2:Nₚ+1,2]
    I_a = sum(I_a_vec)
    I_b = sum(I_b_vec)   
    
    temp_a = S_a*γ_a.(strategies).*(I_a_vec*K_aa/N_a.+I_b_vec*K_ab/N_b)
    temp_b = S_b*γ_b.(strategies).*(I_a_vec*K_ab/N_a.+I_b_vec*K_bb/N_b)
        
    du[1,1] = c*I_a-sum(temp_a)  
    du[1,2] = c*I_b-sum(temp_b)

    du[2:Nₚ+1,1] =temp_a.-c*I_a_vec
    du[2:Nₚ+1,2] =temp_b.-c*I_b_vec 
end

function draw_trait_distribution(sol,timeline,strategies,Nₚ,ylim,img_name)
    
    trait_distribution = @animate for t in timeline       
        ode_sol = sol(t)

        I_a_vec = ode_sol[2:Nₚ+1,1]
        I_b_vec = ode_sol[2:Nₚ+1,2]
        
        plot(strategies,I_a_vec,ylims=(0,ylim))
        plot!(strategies,I_b_vec)
    end
    gif(trait_distribution, "./fig/ex3/"*img_name*".gif", fps = 15) 
end

function run_ex3()
    Nₚ = 200
    z_vec = range(0.0,stop=1.0,length=Nₚ+1)    

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.6
    μ_b,σ²_b,amplitude_b = 0.75,0.0025,0.6 #0.3

    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)
    γ_b(z) = γ_fun(z,μ_b,σ²_b,amplitude_b)    

    #Plot strategy
    plot(z_vec,γ_a.(z_vec),label="Species A",
    xlabel="Strategy",
    ylabel="γ")
    plot!(z_vec,γ_b.(z_vec),label="Species B")

    
    c = 0.1 #Recovery rate

    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3   
    
    N_a = 10^5
    S_a = 0.999*N_a
    I_a = 0.001*N_a
    N_b = 10^5
    S_b = N_b
    I_b = 0.0
    
    I_a_pdf(z) = normal_d(z,μ_a,0.0005)
    I_a_pmf = cont2disc(I_a_pdf,z_vec)    
    
    t_start = 0
    t_end = 600
    tspan = (t_start,t_end)

    u₀ = zeros(Nₚ+1,2)
    u₀[1,1],u₀[1,2]=S_a,S_b 
    u₀[2:Nₚ+1,1] = I_a*I_a_pmf        
 
    p = (c,γ_a,γ_b,z_vec[1:Nₚ],Nₚ,K_aa,K_bb,K_ab,N_a,N_b)
    
    prob = ODEProblem(epievodyn_simple!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE
    
    plot(z_vec[1:Nₚ],sol(150)[2:Nₚ+1,1])
    plot!(z_vec[1:Nₚ],sol(150)[2:Nₚ+1,2])

    timeline = range(t_start,stop=t_end,length=200)

    draw_trait_distribution(sol,timeline,z_vec[1:Nₚ],Nₚ,0.7*N_a,"ex3")

    S_a_fun(t) = sum(sol(t)[1,1])
    S_b_fun(t) = sum(sol(t)[1,2])
    I_a_fun(t) = sum(sol(t)[2:Nₚ+1,1])
    I_b_fun(t) = sum(sol(t)[2:Nₚ+1,2])

    plot(timeline,S_a_fun.(timeline),label="S_a")
    plot!(timeline,S_b_fun.(timeline),label="S_b")
    plot!(timeline,I_a_fun.(timeline),label="I_a")
    plot!(timeline,I_b_fun.(timeline),label="I_b")
end

run_ex3()