#Single-strain-system R₀
function R₀_fun(Δz,Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)
    a  = exp(-Δz^2/(2*σ²))
    b =  exp(-(Δx^2-2*Δz*Δx)/(2*σ²))
    R₀_mean = (R₀_aa_max+R₀_bb_max*b)/2
    R₀_diff = (R₀_aa_max-R₀_bb_max*b)/2

    return a*(R₀_mean+sqrt(R₀_diff^2+R₀_ab_max^2*b))
end 

function calc_diff(x,f;h=0.00001)
    return (f(x+h)-f(x))/h
end

function calc_sec_diff(x,f;h=0.00001)
    return (f(x+h)-2*f(x)+f(x-h))/h^2
end

function find_min_system_R₀(Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²;n_points=100)
    f(x) = R₀_fun(x,Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)   

    Δz_rel = collect(range(0.0,stop=1.0,length=n_points))
    Δz = Δz_rel*Δx

    return minimum(f.(Δz))
end

function run_test()
    c = 0.1 #Intraspecific transmission rate coefficint

    #maximum intraspecific basic reproduction number
    R₀_aa_max = 1.5 
    R₀_bb_max = 1.8 
    R₀_ab_max = (R₀_aa_max+R₀_bb_max)/2*c    
    
    σ² = 0.0025   

    Δx = 0.09

    Δz_rel = collect(range(0.0,stop=1.0,length=100))
    Δz = Δz_rel*Δx

    R₀ = R₀_fun.(Δz,Ref(Δx),Ref(R₀_aa_max),Ref(R₀_bb_max),Ref(R₀_ab_max),Ref(σ²))
    
    plot(Δz,R₀)    

    f(x) = R₀_fun(x,Δx,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)        
    df(x) = calc_diff(x,f)  
    ddf(x) = calc_sec_diff(x,f)      

    ϵ = 0.0001
   
    @show a = 0.0-df(0.0)/ddf(0.0)+ϵ
    @show df(a)
    
    @show b = Δx-df(Δx)/ddf(Δx)-ϵ
    @show df(b)
    
    @show df(Δx-0.00001)

    #@show bisection(df, a, b)
    #plot(Δz,df.(Δz)) 
    #plot!([minimum(Δz),maximum(Δz)],[0.0,0.0])  
    #plot!([a],[df(a)],seriestype=:scatter)  
    #plot!([b],[df(b)],seriestype=:scatter) 
    
    #plot(Δz,f.(Δz)) 

    @show minimum(f.(Δz))

    R₀_min(x) = find_min_system_R₀(x,R₀_aa_max,R₀_bb_max,R₀_ab_max,σ²)

    Δx_vec = collect(range(0.01,stop=0.2,length=100))

    plot(Δx_vec,R₀_min.(Δx_vec))
end

run_test()