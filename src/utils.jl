
function molar_to_weight(value,x,MW)
    return value/LinearAlgebra.dot(MW,x)*1000.0
end

function weight_to_molar(value,x,MW)
    return (value*LinearAlgebra.dot(MW,x))/1000.0
end


#transformation of arbitrary molar units to the SI units
@inline function _transform_v(v,mw,x)
    return float(v)
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐍^-1},mw,x)
    return float(Unitful.ustrip(u"m^3/mol"(v)))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐍^1},mw,x)
    return float(inv(Unitful.ustrip(u"mol/m^3"(v))))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐌^1},mw,x)
    return weight_to_molar(1/Unitful.ustrip(u"kg/m^3"(v)),mw,x)
end


@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐌^-1},mw,x)
    return weight_to_molar(Unitful.ustrip(u"m^3/kg"(v)),mw,x)
end

@inline function _transform_T(T)
    return float(T)
end

@inline function _transform_T(T::TT) where TT <: Unitful.Temperature
        return float(Unitful.ustrip(u"K"(T)))
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

function solve_cubic_eq(coeff1::T,coeff2::T,coeff3::T,coeff4::T) where T<:Real
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    
    a1  =  1 / coeff4
    E1  = -coeff3*a1
    E2  =  coeff2*a1
    E3  = -coeff1*a1
    s0  =  E1
    E12 =  E1*E1
    A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B   =  E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt(A*A - 4*B*B*B)
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2)
end

