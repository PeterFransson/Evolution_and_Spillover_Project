function within_host_pip()
    μ_a,σ²_a,amplitude_a = 0.4,0.0025,0.6
    γ_a(z) = γ_fun(z,μ_a,σ²_a,amplitude_a)

    Nₚ=200
    z_vec = range(0.0,stop=1.0,length=Nₚ)

    
    r_m_mat = zeros(Nₚ,Nₚ)
    for col_idx in 1:1:Nₚ #Resident strain
        resident_strain = z_vec[col_idx]
        for i in 1:1:Nₚ #Mutant strain strain
            mutant_strain = z_vec[i]
            if γ_a(mutant_strain)>=γ_a(resident_strain)
                r_m_mat[i,col_idx]=1
            end            
        end
    end   
    
    heatmap(z_vec,z_vec,r_m_mat)
    plot!([0.0,1.0],[0.0,1.0],legend=false)   
    savefig("./fig/AD_example_1/within_host_PIP.svg") 
end

within_host_pip()