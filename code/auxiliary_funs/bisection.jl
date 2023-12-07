#Bisection method for finding one root to f in the interval [a,b]
function bisection(f::Function, a::AbstractFloat, b::AbstractFloat;
    tol::AbstractFloat=1e-6, maxiter::Integer=200)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
    i += 1
    i != maxiter || error("Max iteration exceeded")
    c = (a+b)/2
    fc = f(c)
    if fc == 0
    break
    elseif fa*fc > 0
    a = c  # Root is in the right half of [a,b].
    fa = fc
    else
    b = c  # Root is in the left half of [a,b].
    end
    end
    return c
end