normal_d(x,μ,σ²) = exp(-(x-μ)^2/(σ²*2))/sqrt(2*π*σ²)

function draw_trait_plot()
    trait_line = range(0.0,stop=1.0,length=200)

    μ_A = 0.2
    μ_B = 0.3
    μ_C = 0.75
    σ² = 0.0025
    scale_fac = sqrt(2*π*σ²)*0.6
    γ_fun(x,μ,σ²) = scale_fac*normal_d(x,μ,σ²)
    plot(trait_line,γ_fun.(trait_line,Ref(μ_A),Ref(σ²)),xlabel="Pathogen Trait",ylabel="γ",label="A")
    plot!(trait_line,γ_fun.(trait_line,Ref(μ_B),Ref(σ²)),label="B")
    plot!(trait_line,γ_fun.(trait_line,Ref(μ_C),Ref(σ²)),label="C")
    savefig("./fig/Trait_Plot/trait_plot.svg")
end

draw_trait_plot()