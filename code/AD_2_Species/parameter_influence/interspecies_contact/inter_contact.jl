#Run computation
function compute() 
    δz_vec = range(0.05,stop=0.15,length=50)
    c_vec = [0.2 0.3 0.4 0.5 0.6 0.7 0.8] #[0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8]

    γ = 0.1 
    z_start = 0.18
    σ²,amplitude = 0.0025,0.6
    μ_a = 0.2
    Nₚ = 1500    

    #Parameters
    c_aa = 0.425
    c_bb = 0.425
    
    
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
    Threads.@threads for i in eachindex(c_vec)
        for j in eachindex(δz_vec)
            δz =  δz_vec[j]
            c = c_vec[i]
            file_name = "output_c_$(i)_dz_$(j)"
            
            c_ab = (c_aa+c_bb)/2*c
            μ_b = 0.2+δz        
            z_end = 0.2+δz+0.05

            option = z_start,z_end,Nₚ,u₀,tspan

            syspar = SystemParameters(μ_a,μ_b,c_aa,c_bb,c_ab,N_a,N_b,γ,γ,σ²,amplitude)

            strats = singular_strategies(syspar,option)
            
            JLD2.@save "./output/AD_2_Species/parameter_influence/interspecies_contact/"*file_name*".jld2" strats syspar
        end
    end

    c_len = length(c_vec)
    dz_len = length(δz_vec)

    JLD2.@save "./output/AD_2_Species/parameter_influence/interspecies_contact/c_dz_info.jld2" c_len dz_len
end

function get_n_types(strats::Vector{SingularStrat})
    if isempty(strats)
        return nothing
    else
        ntypes = zeros(Integer,4)
        for strat in strats
            if strat.conv_stable&&strat.evo_stable
                ntypes[1] += 1
            elseif strat.conv_stable&&!strat.evo_stable
                ntypes[2] += 1
            elseif !strat.conv_stable&&strat.evo_stable
                ntypes[3] += 1
            else
                ntypes[4] += 1
            end
        end
        return ntypes
    end
end

function find_type(strats::Vector{SingularStrat})

    ntypes = get_n_types(strats)

    if length(strats)==1&&ntypes[1]==1
        return :blue       
    elseif length(strats)==1&&ntypes[2]==1
        return :green
    elseif length(strats)==3&&ntypes[1]==2&&ntypes[4]==1
        return :red
    else
        return :black
    end
end

type_dict = Dict(:blue =>(1,"Type A"),:green=>(2,"Type B"),:red=>(4,"Type C"),:black=>(3,"Transient Type"))

function draw_plot()
    c_len = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/c_dz_info.jld2","c_len")
    dz_len = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/c_dz_info.jld2","dz_len")
    fig_name = "interspecies_contact"

    type_width = ((Float64[],Float64[],:blue),
    (Float64[],Float64[],:green),
    (Float64[],Float64[],:black),
    (Float64[],Float64[],:red))

    plot(xlabel="Between Species Distance",ylabel="Interspecies Contact")

    for i = 1:c_len
        file_name = "output_c_$(i)_dz_$(1)" 
        syspar = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/"*file_name*".jld2","syspar")
        strats = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/"*file_name*".jld2","strats")
        c_ab = syspar.c_ab 
        c_aa = syspar.c_aa
        c_bb = syspar.c_bb
        c = c_ab*2/(c_aa+c_bb)
        δz = abs(syspar.z_a-syspar.z_b)
        type_sym = find_type(strats)
        for j = 2:dz_len 
            file_name_temp = "output_c_$(i)_dz_$(j)" 
            strats_temp = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/"*file_name_temp*".jld2","strats")
            syspar_temp = JLD2.load("./output/AD_2_Species/parameter_influence/interspecies_contact/"*file_name_temp*".jld2","syspar")
            δz_temp = abs(syspar_temp.z_a-syspar_temp.z_b)
            type_sym_temp = find_type(strats_temp)
            if (type_sym != type_sym_temp)||j==dz_len
                if i == 1
                    plot!([δz,δz_temp],[c,c],color=type_sym,markershape=:circle,label=type_dict[type_sym][2])
                else
                    plot!([δz,δz_temp],[c,c],color=type_sym,markershape=:circle,label="")
                end

                push!(type_width[type_dict[type_sym][1]][1],c)
                push!(type_width[type_dict[type_sym][1]][2],abs(δz-δz_temp))
                
                δz = δz_temp    
                type_sym = type_sym_temp          
            end
        end
    end   
    plot!(legend = false) 
       
    savefig("./fig/AD_2_Species/parameter_influence/interspecies_contact/"*fig_name*"_subpl_1.svg")

    plot(xlabel="Interspecies Contact",ylabel="Type Width")
    for yx in type_width
        plot!(yx[1],yx[2],color=yx[3],markershape=:circle,label=type_dict[yx[3]][2])
    end
    plot!(legend = :outerbottomright)
    
    savefig("./fig/AD_2_Species/parameter_influence/interspecies_contact/"*fig_name*"_subpl_2.svg")
end
compute()
draw_plot()
#Draw gif
