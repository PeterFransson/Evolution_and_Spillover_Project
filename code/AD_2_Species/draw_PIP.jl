function plot_list()
    folder_name = "./fig/PIP/2_species/pop_test"
    name = "test_3"
    isdir(folder_name)|| mkdir(folder_name)

    Nₚ = 1000
    z_start = 0.1
    z_end = 0.5

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.35
    
    N_tot = 2*10^3 
    r_N = 0.1
    N_a = N_tot*r_N#10^3  
    N_b = N_tot*(1-r_N)#10^3

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 5500
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]    

    option = z_start,z_end,Nₚ,u₀,tspan 

    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,μ_a,μ_b
    
    fig_name = folder_name*"/"*name    
    draw_PIP(fig_name,p_in,option)
    p_system = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)
    savesystemparameters(p_system,folder_name*"/"*name*"_param")
end
plot_list()