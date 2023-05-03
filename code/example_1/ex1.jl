#Simple SI-model example with one species and one adaptive pathogen trait (β)
#Here we asume the trait and mutation distribution are normal distributed 

#ODE system describing the dynamics of the epi-evolutionary model 
#u[1] = S, Susceptible individuals
#u[2] = I, Infected individuals
#u[3] = Eβ, Mean infection rate, β, in the infected population, β~N(Eβ,σ²β)
#u[4] = σ²β, Variance of β~N(Eβ,σ²β)
function epievodyn_simple!(du,u,p,t)
    γ,N₀,μₘ,δₘ,σ²ₘ = p
    S,I,Eβ,σ²β = u
    du[1] = I*γ-S*Eβ*I/N₀
    du[2] = -du[1]
    du[3] = S/N₀*σ²β+μₘ*δₘ   
    du[4] = μₘ*(δₘ^2+σ²ₘ)
end

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)

function draw_trait_distribution(trait_μ,trait_σ,trait_line,img_name)
    length(trait_μ)==length(trait_σ)||error("length(trait_μ)≠length(trait_σ)")
    σ²_min = minimum(trait_σ)
    ymax = 1.01/(sqrt(2*π*σ²_min)) 
    trait_distribution = @animate for (μ,σ²) in zip(trait_μ,trait_σ)        
        plot(trait_line,normal_d.(trait_line,Ref(μ),Ref(σ²)),ylims=(0.0,ymax))
    end
    gif(trait_distribution, "./fig/ex1/"*img_name*".gif", fps = 15) 
end

function run_ex1()
    γ = 0.3 #Recovery rate    
    Eβ₀ = 0.3 #Mean infection rate, β, in the infected population, β~N(Eβ,σ²β)
    σβ₀  = 0.05 #Standard deviation of β~N(Eβ,σ²β) 
    σ²β₀  = σβ₀^2 #Variance of β~N(Eβ,σ²β)
    σ²ₘ = 0.01^2 #Variance of the mutan trait distribution within strain i infected individual
                 #βᵢₘ~N(βᵢ+δₘ,σ²ₘ)
    δₘ = 0.01 #The offset from  he mutan trait distribution within strain i infected individual
              #βᵢₘ~N(βᵢ+δₘ,σ²ₘ)
    μₘ = 0.05 #mutation rate   

    N₀ = 10^5 #Population size
    f_I = 0.001 #Initial fraction of infected individuals
    S₀ = N₀*(1-f_I) #Initial number of susceptible individuals
    I₀ = N₀*f_I #Initial number of infected individuals

    t_end = 100
    tspan = (0.0,t_end) #Solve ODE from 0.0 to t_end
    u₀ = [S₀,I₀,Eβ₀,σ²β₀] #Inital state for the ODE
    p = (γ,N₀,μₘ,δₘ,σ²ₘ) #Paramters for the ODE

    prob = ODEProblem(epievodyn_simple!,u₀,tspan,p) #Setup the ODE problem
    sol = solve(prob) #Solve the ODE problem, sol contains a continuous approximation to the ODE

    get_result(tvec,idx) = [sol(t)[idx] for t in tvec]   
    
    t_line = range(0.0,stop=t_end,length=200)
    
    β_vec = range(0.0,stop=1.0,length=200) #Points along the trait axis of β, used for ploting
                                           #the distribution of β 
    Eβ_vec = get_result(t_line,3)
    σ²β_vec = get_result(t_line,4)

    draw_trait_distribution(Eβ_vec,σ²β_vec,β_vec,"ex1")    
        
    plot(t_line,get_result(t_line,1)/N₀,label="S")
    plot!(t_line,get_result(t_line,2)/N₀,label="I")    
    plot!(t_line,get_result(t_line,3),label="Eβ")
end

run_ex1()