
function molar_to_weight(value,x,MW)
    return value/LinearAlgebra.dot(MW,x)*1000.0
end

function weight_to_molar(value,x,MW)
    return (value*LinearAlgebra.dot(MW,x))/1000.0
end


#transformation of arbitrary molar units to the SI units
@inline function _transform_v(v,mw,x)
    return 1.0*v
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐍^-1},mw,x)
    return Unitful.ustrip(u"m^3/mol"(v))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐍^1},mw,x)
    return 1/Unitful.ustrip(u"mol/m^3"(v))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐌^1},mw,x)
    return weight_to_molar(1/Unitful.ustrip(u"kg/m^3"(v)),mw,x)
end


@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐌^-1},mw,x)
    return weight_to_molar(Unitful.ustrip(u"m^3/kg"(v)),mw,x)
end

@inline function _transform_T(T)
    return 1.0*T
end

@inline function _transform_T(T::TT) where TT <: Unitful.Temperature
        return Unitful.ustrip(u"K"(T))
end

function _transformVT(V,T,mw,x)
    return (_transform_v(V,mw,x),_transform_T(T))
end

@inline function _transform_P(P)
    return 1.0*P
end

@inline function _transform_P(P::TT) where TT <: Unitful.Pressure
        return Unitful.ustrip(u"Pa"(P))
end


###
### fraction utilities
###
function normalizefrac!(x)
    summ = sum(x)
    for i = 1:length(x)
        x[i] /=summ
    end
    x
end

function normalizefrac(x)
    summ = sum(x)
    y = copy(x)
    for i = 1:length(x)
        y[i] /=summ
    end
    return y
end
function randfrac(N::Integer)
    x = rand(N)
    return x ./sum(x)
end
function randfrac(x0::Array)
    length(x0)>0 && begin
    x = similar(x0)
    zerox = 0.0
    for i = 1:length(x0)
     x[i]=ifelse(x0[i]==zerox,zerox,rand())
    end
    x /= sum(x)
    return x 
end
end

function randpurefrac(x0::Array)
    x = zeros(length(x0))
    i = rand(LinearIndices(findall(x->x!=0.0,x0)))
    x[i] = 1.0
    return x
end

