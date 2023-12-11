struct SingularStrat
    strategy::Real #Position in compability space of evolutionarily singular strategy
    conv_stable::Bool #Is the strategy convergence stable strategy
    evo_stable::Bool #Is the strategy evolutionarily stable strategy
end

function calc_selectgrad(strain::Real,p_in,option)
    u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b = p_in

    p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
    eq_point = find_eq(u₀,tspan,p_eq)

    p_grad = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain)
    S = (eq_point.S_a,eq_point.S_b)
    dλ_max = calculate_select_grad(p_grad,S)

    return dλ_max
end
function calc_d_selectgrad(strain::Real,p_in,option)
    u₀,tspan = option    
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b = p_in

    p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
    eq_point = find_eq(u₀,tspan,p_eq)

    p_curv = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain),τ_d_prime_a(strain),τ_d_prime_b(strain)
    S = (eq_point.S_a,eq_point.S_b)
    ddλ_max = calculate_sec_inv_fit(p_curv,S)

    return ddλ_max
end

function check_convergence(strain::Real,p_in,option;h::Real=10^-3) where {R<:Real,T<:Real}
    
    dλ_max_1 = calc_selectgrad(strain+h,p_in,option)
    dλ_max_2 = calc_selectgrad(strain-h,p_in,option)

    return (dλ_max_1-dλ_max_2)/(2*h)
end

function calc_selectgrad(p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a,τ_prime_b = p_in

    strategies = range(z_start,stop=z_end,length=Nₚ) 
    select_grad = zeros(Nₚ)

    for i = 1:Nₚ
        strain = strategies[i]
        p_eq = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b 
        eq_point = find_eq(u₀,tspan,p_eq)

        p_grad = τ_a(strain),τ_b(strain),c_aa,c_bb,c_ab,γ_a,γ_b,N_a,N_b,τ_prime_a(strain),τ_prime_b(strain)
        S = (eq_point.S_a,eq_point.S_b)
        dλ_max = calculate_select_grad(p_grad,S)

        select_grad[i] = dλ_max
    end

    return (select_grad,collect(strategies))
end

function find_bracketing_interval(yvals::Vector{R},xvals::Vector{T}) where {R<:Real,T<:Real}
    N = length(yvals)
    y_signs = yvals.<0.0    
    if sum(y_signs)==0||sum(y_signs)==N
        println("No sign change detected in the given intervall")
        return nothing        
    else
        bracketing = Tuple{Real,Real}[]
        left_sign = y_signs[1]
        left_val = xvals[1]
        for i = 2:N 
            right_sign = y_signs[i]
            if right_sign!=left_sign
                right_val = xvals[i]
                push!(bracketing,(left_val,right_val)) 

                left_sign = right_sign
                left_val = right_val
            end
        end
        return bracketing
    end
end
    
function find_singular_strategies(p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option    
    option_bisec = u₀,tspan
    Nₚ>2||error("Nₚ≤2")
    select_grad,strategies = calc_selectgrad(p_in,option)    
    b_vals = find_bracketing_interval(select_grad,strategies)
    isnothing(b_vals)==false||error("No sign change detected in the given intervall")
    z_0s = zeros(length(b_vals))    
    for (i,val) in enumerate(b_vals)        
        z_0s[i] = bisection(x->calc_selectgrad(x,p_in,option_bisec), val[1],val[2])             
    end
    return z_0s
end

function singular_strategies(p_in,option)
    z_start,z_end,Nₚ,u₀,tspan = option
    τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b = p_in

    p = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b
    option_in = u₀,tspan

    Δz = (z_end-z_start)/(Nₚ-1)

    z0s = find_singular_strategies(p,option)     

    sing_strats = SingularStrat[]

    for z0 in z0s         
        conv_check = check_convergence(z0,p,option_in;h=Δz)    
        evo_check = calc_d_selectgrad(z0,p_in,option_in)
        push!(sing_strats,SingularStrat(z0,conv_check<0,evo_check<0))
    end

    return sing_strats
end

function plot_def()
    Nₚ = 1500
    z_start = 0.18
    z_end = 0.33
    Δz = (z_end-z_start)/(Nₚ-1)

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    c_ab = 0.3 #0.3

    γ = 0.1 
    σ²,amplitude = 0.0025,0.6
    μ_a,μ_b = 0.2,0.323
    
    N_a = 10^3   
    N_b = 10^3 

    #Initial values    
    I_a₀ = 1  
    S₀_a = N_a-I_a₀ 
         
    I_b₀ = 1
    S₀_b = N_b-I_b₀     

    τ_a(z) = τ_fun(z,μ_a,σ²,amplitude)
    τ_b(z) = τ_fun(z,μ_b,σ²,amplitude)
    τ_prime_a(z) = τ_prime_fun(z,μ_a,σ²,amplitude)
    τ_prime_b(z) = τ_prime_fun(z,μ_b,σ²,amplitude)
    τ_d_prime_a(z) = τ_d_prime_fun(z,μ_a,σ²,amplitude)
    τ_d_prime_b(z) = τ_d_prime_fun(z,μ_b,σ²,amplitude)

    t_start = 0
    t_end = 7000
    tspan = (t_start,t_end)

    u₀ = [S₀_a I_a₀;S₀_b I_b₀]     

    option = z_start,z_end,Nₚ,u₀,tspan
    p_in = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b
    p_strats = τ_a,τ_b,c_aa,c_bb,c_ab,γ,γ,N_a,N_b,τ_prime_a,τ_prime_b,τ_d_prime_a,τ_d_prime_b

    select_grad,strategies = calc_selectgrad(p_in,option)

    strats = singular_strategies(p_strats,option)
    @show strats 

    #=
    option_bisec = u₀,tspan
    
    z0s = find_singular_strategies(p_in,option)
    @show z0s
    conv_vecs = check_convergence.(z0s,Ref(p_in),Ref(option_bisec);h=Δz)
    @show conv_vecs
    ddλ_max_0s = calc_d_selectgrad.(z0s,Ref(p_d_selectgrad),Ref(option_bisec))
    @show ddλ_max_0s
    =#

    plot([z_start,z_end],[0,0],ylim=[-0.1,0.1])
    plot!(strategies,select_grad,ylabel="Sᵣ(m=r)′",xlabel="Strain",legends=false)
    #savefig("./fig/PIP/2_species/between_case_b_c/singular_strategies.svg")
end

plot_def()