module Taylor9


export Taylor,paso,horner,serieTaylor,taylorIntegration,Picard

mutable struct Taylor{T<:Real}
    orden::Int64
    coef::Array{T, 1}
    Taylor(coef::Array{T, 1}) = new(length(coef)-1, coef)
end

import Base: +, -, *, /,^, getindex, setindex!,==

function getindex(f::Taylor, i::Int64)
    f.coef[i + 1]
end

function getindex(f::Taylor, r::UnitRange{Int64})
    newcoef = f.coef[(r[1] +1): (r[end] +1)]
    Taylor(newcoef)
end

function setindex!(f::Taylor, x::Float64, i::Int64)
    f.coef[i + 1] = x
end

function setindex!(f::Taylor, x::Int64, i::Int64)
    f.coef[i + 1] = x
end

function ==(f::Taylor, g::Taylor)
    f.orden == g.orden && f.coef == g.coef  ? true : false
end

function +(f::Taylor, g::Taylor)
    if f.orden == g.orden
        Taylor(f.coef .+ g.coef)
    elseif f.orden < g.orden
        Taylor(f.coef[1:f.orden + 1] .+ g.coef[1:f.orden + 1])
    else
        Taylor(g.coef[1:g.orden + 1] .+ f.coef[1:g.orden + 1])
    end
end

function -(f::Taylor, g::Taylor)
    if f.orden == g.orden
        Taylor(f.coef .- g.coef)
    elseif f.orden < g.orden
        Taylor(f.coef[1:f.orden + 1] .- g.coef[1:f.orden + 1])
    else
        Taylor(-g.coef[1:g.orden + 1] .+ f.coef[1:g.orden + 1])
    end
end

function *(f::Taylor, g::Taylor)
    f.orden != g.orden && error("Ordenes distintos")
    newcoef = Array{Float64, 1}(f.orden + 1)
    for k in 0:f.orden 
        Taylori = 0.0
        for i in 0:k
            Taylori += f[i]*g[k - i]  
        end
        newcoef[k + 1] = Taylori 
    end
    Taylor(newcoef)
end

function /(f::Taylor, g::Taylor)
    f.orden != g.orden && error("Ordenes distintos")
    newcoef = Array{Float64, 1}(f.orden + 1)
    
    for k in 0:f.orden 
        Taylori = f[k]
        for i in 0:k - 1
            Taylori -= (f[i]/g[i])*g[k - i]  
        end
        newcoef[k + 1] = Taylori/g[0] 
    end
    Taylor(newcoef)
end

function ^(f::Taylor, n::Int64)
    newTaylor = f
    for i in 1:n-1
        newTaylor *= f
    end
    newTaylor
end

function paso(x::Taylor)
    p = length(x.coef) - 1
    (1.0e-32/(abs(x[p])))^(1/p)
end

function horner(P::Taylor, x::Real)
    result = 0.0
    for i in P.orden:-1:0
        result = (result * x) + P[i]
    end
    result
end

function serieTaylor(f,  x0::Real, t0::Real, orden::Int64=10)
    newTaylor = Taylor(zeros(orden + 1))
    newTaylor[0] = x0
    newTaylor[1] = f(x0, t0)
    for i in 2:orden
        x0 = newTaylor[0:i]
        newTaylor[i] = (f(x0, t0)[i-1])/(i)
    end
    newTaylor
end

function taylorIntegration(f, x0::Real, t0::Real ,orden::Int64=32)
    T = serieTaylor2(f, x0, t0, orden)
    t = paso(T)
    x = horner(T, t)
    x, t + t0
end

function Picard(f, x0::Real, t::Real, orden::Int64=32)
    x0new = Taylor(zeros(orden + 1))
    x0new[0] = x0
    x0new[1] = f(x0)
    for i in 2: orden 
        x0 = x0new[0:i]
        x0new[i] = (f(x0)[i-1])/(i)
    end
    horner(x0new, t)
end

end