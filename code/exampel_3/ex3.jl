#Simple SIS-model example with two species and one strategy (z:=[0,1])

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
    I_a_vec = u[1,1:Nₚ]
    I_b_vec = u[2,1:Nₚ]
    I_a = sum(I_a_vec)
    I_b = sum(I_b_vec)
    
    f_a = I_a_vec.*γ_a.(strategies)*K_aa/N_a
    f_b = I_b_vec.*γ_b.(strategies)*K_aa/N_a
        
    du[1,1] = c*I_a-S_a*sum()      
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

function run_ex3()
    z_vec = range(0.0,stop=1.0,length=300)

    μ_a,σ²_a,amplitude_a = 0.2,0.0025,0.4
    μ_b,σ²_b,amplitude_b = 0.3,0.0025,0.4

    c = 0.1 #Recovery rate

    #Plot strategy
    plot(z_vec,γ_fun.(z_vec,Ref(0.2),Ref(0.0025),Ref(0.6)),label="Species A",
    xlabel="Strategy",
    ylabel="γ")
    plot!(z_vec,γ_fun.(z_vec,Ref(0.3),Ref(0.0025),Ref(0.6)),label="Species B")

    K_aa = 0.425
    K_bb = 0.425
    K_ab = 0.3   
    
    N_a = 10^5
    S_a = 0.999*N_a
    I_a = 0.001*N_a
    N_b = 10^5
    S_b = N_b
    I_b = 0.0
end

run_ex3()