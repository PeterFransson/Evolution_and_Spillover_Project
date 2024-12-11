function get_depressed_cubic_coef(a,b,c,d)
    p = (3*a*c-b^2)/(3*a^2)
    q = (2*b^3-9*a*b*c+27*a^2*d)/(27*a^3)
    return (p,q)
end

#Number used to calcualted the real solution to a depressed cubic equation when calc_discriminant_depressed_cubic>0  
function calc_depressed_cubic_num(p,q)
    return q^2/4+p^3/27
end

#Cardon's formula: 
#Finds the (only) real root of the the Depressed cubic t^3+pt+q=0 
#where discriminant_test(p,q) is true 
function cubic_sol_one_real(p,q)
    depressed_cubic_num = calc_depressed_cubic_num(p,q)

    u₁ = -q/2+sqrt(depressed_cubic_num)
    u₂ = -q/2-sqrt(depressed_cubic_num)

    return cbrt(u₁)+cbrt(u₂)
end

#Finds the three real distinct root of the the Depressed cubic t^3+pt+q=0 
#where discriminant_test(p,q) is false (François Viète)
function cubic_sol_three_real(p,q)
    t₁ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3)
    t₂ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3-2*π*1/3)
    t₃ = 2*sqrt(-p/3)*cos(acos(3*q*sqrt(-3/p)/(2*p))/3-2*π*2/3)

    return (t₁,t₂,t₃)
end
#Finds all distinct real roots of the depressed cubic t^3+pt+q=0
function cubic_sol(p,q)
    depressed_cubic_num = calc_depressed_cubic_num(p,q)
    if depressed_cubic_num≈0
        if p≈0 #One Triple root
            return 0
        else #One double root
            t = -3*q/(2*p)
            return (3*q/p,t)            
        end    
    elseif depressed_cubic_num>0
        return cubic_sol_one_real(p,q)
    else
        return cubic_sol_three_real(p,q)
    end
end

#Finds all distinct real solutions to ax³+bx²+cx+d=0 when (a,b,c,d) are real
function cubic_sol(a,b,c,d)
    (p,q) = get_depressed_cubic_coef(a,b,c,d)
    
    t_sols = cubic_sol(p,q)
    return t_sols.-b/(3*a)
end