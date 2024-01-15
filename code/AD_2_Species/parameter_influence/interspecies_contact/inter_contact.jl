#Run computation
function compute() 
    δz_vec = range(0.05,stop=0.15,length=20)

    γ = 0.1 
    z_start = 0.18
    σ²,amplitude = 0.0025,0.6
    μ_a = 0.2
    Nₚ = 1500    

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3
    
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]      
        
    μ_b = 0.2+δz        
    z_end = 0.2+δz+0.05

    option = z_start,z_end,Nₚ,u₀,tspan
    syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)
end

#Draw gif
