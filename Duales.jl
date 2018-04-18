module Duales

export Dual, derivada, grad

struct Dual{T<:Real}
    fa  :: T
    fpa :: T
end

import Base: +, -, *, /, exp 

function +(x::Dual, y::Dual)
    Dual(x.fa + y.fa, x.fpa + y.fpa)
end

function +(c::Real, y::Dual)
    Dual(c + y.fa, y.fpa)
end

function +(x::Dual, c::Real)
    Dual(x.fa + c, x.fpa )
end


function -(x::Dual, y::Dual)
    Dual(x.fa - y.fa, x.fpa - y.fpa)
end

function -(c::Real, y::Dual)
    Dual(c - y.fa, y.fpa)
end

function -(x::Dual, c::Real)
    Dual(x.fa - c, x.fpa )
end

function *(x::Dual, y::Dual)
    Dual(x.fa * y.fa, x.fpa * y.fa + x.fa * y.fpa )
end

function *(c::Real, y::Dual)
    Dual(c* y.fa,  c * y.fpa )
end

function *(x::Dual, c::Real)
    Dual(x.fa*c, x.fpa*c )
end

function ^(x::Dual, n::Int64)
    Dual(x.fa^n, n*x.fa^(n-1))
end

function exp(x::Dual)
    Dual(exp(x.fa), exp(x.fa)*x.fpa)
end

function derivada(f, x0::Real)
    x = f(Dual(x0, one(x0)))
    x.fpa
end

function grad(f,x0,y0)
    xaux = f(Dual(x0, 1), y0)
    yaux = f(x0, Dual(y0, 1))

    Df = (xaux.fpa,yaux.fpa)

end

end