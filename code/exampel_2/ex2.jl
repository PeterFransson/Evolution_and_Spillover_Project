#Simple SIS-model example with two species and one adaptive pathogen trait (z:=[0,1])
#Here we asume the trait and mutation distribution are normal distributed 

normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)
γ_fun(x,μ,σ²,amplitude) = amplitude*exp(-(x-μ)^2/(σ²*2))

function run_ex2()
    trait_line = range(0.0,stop=1.0,length=200)
    μ_A = 0.2
    μ_B = 0.3
    σ² = 0.0025

    γ = 0.3 #Recovery rate    
    Eβ₀ = 0.3 #Mean infection rate, β, in the infected population, β~N(Eβ,σ²β)
    σβ₀  = 0.05 #Standard deviation of β~N(Eβ,σ²β) 
    σ²β₀  = σβ₀^2 #Variance of β~N(Eβ,σ²β)
    σ²ₘ = 0.01^2 #Variance of the mutan trait distribution within strain i infected individual
                 #βᵢₘ~N(βᵢ+δₘ,σ²ₘ)
    δₘ = 0.01 #The offset from  he mutan trait distribution within strain i infected individual
              #βᵢₘ~N(βᵢ+δₘ,σ²ₘ)
    μₘ = 0.05 #mutation rate   

    plot(trait_line,γ_fun.(trait_line,Ref(μ_A),Ref(σ²),Ref(0.6)),xlabel="Pathogen Trait",ylabel="γ",label="A")
    plot!(trait_line,γ_fun.(trait_line,Ref(μ_B),Ref(σ²),Ref(0.6)),label="B")
end

run_ex2()