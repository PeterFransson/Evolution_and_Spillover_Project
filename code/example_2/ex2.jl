#Simple SIS-model example with one species and one adaptive pathogen trait (z:=[0,1])

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
    c,N₀,μₘ,K_AA,Nₚ,z_vec,μ_a,σ²_a,A_a,Mₘ = p
    S,I = u[1],u[2]
    I_vec = u[3:2+Nₚ]      
    du[1] = I*c-S*sum(K_AA/N₀*γ_fun.(z_vec,Ref(μ_a),Ref(σ²_a),Ref(A_a)).*I_vec) #S
    du[2] = -du[1]     
    for i in 1:Nₚ 
        t1 = S*K_AA/N₀*γ_fun(z_vec[i],μ_a,σ²_a,A_a)*I_vec[i]
        t2 = -c*I_vec[i]
        t3 = -μₘ*I_vec[i]*sum(ρ_fun.(z_vec,Ref(I_vec[i]),Ref(μ_a),Ref(σ²_a),Ref(A_a))*Mₘ)
        t4_prim = [ρ_fun(I_vec[i],I_vec[j],μ_a,σ²_a,A_a)*Mₘ*I_vec[j] for j in 1:Nₚ]
        t4 = μₘ*sum(t4_prim)        
        du[i+2]= t1+t2+t3+t4
    end         
end

function draw_trait_distribution(sol,timeline,trait_line,img_name)
    ymin = 0
    ymax = 1.0
    xmin = minimum(trait_line)
    xmax = maximum(trait_line)
    trait_distribution = @animate for t in timeline       
        ode_sol = sol(t)
        I = ode_sol[2]
        freq = ode_sol[3:end]/I
        plot(trait_line,freq,seriestype=:scatter,ylims=(ymin,ymax),xlims=(xmin,xmax))
    end
    gif(trait_distribution, "./fig/ex2/"*img_name*".gif", fps = 15) 
end

function run_ex2()
    Nₚ = 100
    Δₚ = 1/(Nₚ-1)
    z_target = 0.4
    @show z_idx = ceil(Int64,z_target/Δₚ+1)    

    trait_line = range(0.0,stop=1.0,length=Nₚ)
    μ_A = 0.2    
    σ² = 0.0025
    
    c = 0.3 #Recovery rate     
    μₘ = 0.05 #mutation rate   

    K_AA = 0.7
           
    Mₘ = 1/Nₚ
    N₀ = 10^5
    I_ratio = 0.01
    I = I_ratio*N₀
    S = (1-I_ratio)*N₀
    u₀ = zeros(Nₚ+2)
    u₀[1] = S
    u₀[2] = I 
    u₀[3:2+z_idx] .= I/z_idx

    t_end = 100
    tspan = (0.0,t_end) #Solve ODE from 0.0 to t_end
    p = (c,N₀,μₘ,K_AA,Nₚ,trait_line,μ_A,σ²,0.6,Mₘ)

    prob = ODEProblem(epievodyn_simple!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    get_result(tvec,idx) = [sol(t)[idx] for t in tvec]     

    tvec = range(0.0,stop=t_end,length=200)
    draw_trait_distribution(sol,tvec,trait_line,"freq_evolution")

    plot(trait_line,γ_fun.(trait_line,Ref(μ_A),Ref(σ²),Ref(0.6))*K_AA/c,xlabel="Pathogen Trait",ylabel="R0",label="A")
end

run_ex2()