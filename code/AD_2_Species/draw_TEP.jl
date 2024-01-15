function draw_TEP()
    folder_name = "./fig/AD_2_Species_Alt/TEP/between_case_b_c_2"
    name = "between_case_b_c_2"
    isdir(folder_name)|| mkdir(folder_name)

    Nₚ = 1500
    z_start = 0.1
    z_end = 0.34

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.324    
    
    N_a = 10^3   
    N_b = 10^3 

    syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀       

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    
    u₀_2_strain = [N_a-2 1 1;N_b-2 1 1]    
    
    option = z_start,z_end,Nₚ,u₀,tspan 
    
    coex_region = create_coex_region(syspar,option)   
    z_vec = range(z_start,stop=z_end,length=Nₚ) 
    heatmap(z_vec,z_vec,coex_region)     
    
    option = z_start,z_end,Nₚ,u₀,tspan,u₀_2_strain    
    tep = create_TEP(syspar,option,coex_region;ϵ=1e-2)
    
    #Check cone of invasions  
    @show calc_invasion_cone(0.174,0.348,syspar,option)          
    @show calc_invasion_cone(0.178,0.286,syspar,option)   
    @show calc_invasion_cone(0.229,0.340,syspar,option)    
    @show calc_invasion_cone(0.228,0.289,syspar,option)    
    
    heatmap(z_vec,z_vec,tep) 
    plot!([z_start,z_end],[z_start,z_end],c=:green)
    plot!([μ_a,μ_a],[z_start,z_end],c=:red)
    plot!([z_start,z_end],[μ_a,μ_a],c=:red)
    plot!([z_start,z_end],[μ_b,μ_b],c=:blue)
    plot!([μ_b,μ_b],[z_start,z_end],legend=false,c=:blue)
    savefig(folder_name*"/"*name*".svg")     
end

draw_TEP()