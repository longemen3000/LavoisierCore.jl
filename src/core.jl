using JuliaDB, ForwardDiff, DiffResults, Unitful, LinearAlgebra
include("utils.jl")

abstract type AbstractThermoModel end #base model
abstract type AbstractHelmholtzModel <:  AbstractThermoModel end
abstract type  AbstractGibbsModel <:  AbstractThermoModel end
abstract type  AbstractActivityModel <:  AbstractThermoModel end   
abstract type AbstractMixingRule end # for mixing rules


"""
core_helmholtz(model,V,T,x)

calculates the helmholtz energy. uses the following units:
v:molar volume, m3/mol ,SI
T:Temperature, K , SI
x:molar fraction, adimensional


"""
function core_helmholtz(model::M,V,T,x) where M <:AbstractHelmholtzModel
end

function helmholtz(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return core_helmholtz(model,v2,T2,x)*1.0u"J/mol"
end


#
#Stability test
#should return true if the phase is stable, false if not
#



function _gradientvt(model::AbstractHelmholtzModel,v,T,x) 
    d = [v T]
    f(d) = core_helmholtz(model,d[1],d[2],x)
    return ForwardDiff.gradient(f,d)
end

function _gradientvt2(model::AbstractHelmholtzModel,v,T,x) 
    vt = [v T]
    df = DiffResults.GradientResult(vt)
    df = ForwardDiff.gradient!(df,z->core_helmholtz(model,z[1],z[2],x),vt)
    return (DiffResults.value(df),DiffResults.gradient(df))
end

function _dfdv(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->core_helmholtz(model,z,T,x),v)
end

function _dfdt(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->core_helmholtz(model,v,z,x),T)
end

function _hessianvt(model::AbstractHelmholtzModel,v,T,x)
    f(d) = core_helmholtz(model,d[1],d[2],x)
    return ForwardDiff.hessian(f,[v T])
end

function _hessianvt2(model::AbstractHelmholtzModel,v,T,x)
    vt = [v T]
    fhelmholtz = z->core_helmholtz(model,z[1],z[2],x)
    df = DiffResults.HessianResult(vt)
    df = ForwardDiff.hessian!(df,fhelmholtz,vt)
    return (DiffResults.value(df),
    DiffResults.gradient(df),
    DiffResults.hessian(df))
end



struct HelmholtzResultData
    value::Array{Float64,1}
    diffresult::DiffResults.DiffResult
end

function _diffdata(model::AbstractHelmholtzModel ,V,T,x,order::Int64 = 2)
    d =  vcat(V,T,x[:])
    #f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    res = DiffResults.HessianResult(d)
    if order == 0
    elseif order == 1
    res= ForwardDiff.gradient!(res, f, d)
    elseif order == 2
    res= ForwardDiff.hessian!(res, f, d)
    else
        throw(DomainError(order, """

diffdata can't calculate the derivatives of order $order.
only gradients (order == 1) and hessians (order == 2) are implemented.

"""))
    end
    return HelmholtzResultData(d,res)
end
function diffdata(model::AbstractHelmholtzModel ,v,T,x,order::Int64 = 2)
(v2,T2)=_transformVT(v,T,model.molecularWeight,x)
return _diffdata(model,v2,T2,x,order)
end

_valuevt(df::HelmholtzResultData)  = df.value
_helmholtzvt(df::HelmholtzResultData)  = DiffResults.value(df.diffresult)
_gradientvt(df::HelmholtzResultData)  = DiffResults.gradient(df.diffresult)[1:2]
_hessianvt(df::HelmholtzResultData)  = DiffResults.hessian(df.diffresult)[1:2,1:2]

#
#Pressure, a single derivative is faster
#
function _pressure_impl(model::AbstractHelmholtzModel,v,T,x)
    return -_dfdv(model,v,T,x)
return 
end

function _pressure(df::HelmholtzResultData)
    return -_gradientvt(df)[1]
return 
end

function pressure(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _pressure_impl(model,v2,T2,x)*1.0u"Pa"
end

function pressure(df::HelmholtzResultData) 
    return _pressure(df)*1.0u"Pa"
end

function compressibility_factor(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _pressure_impl(model,v2,T2,x)*v2/(ustrip(Unitful.R)*T2)
end

function compressibility_factor(df::HelmholtzResultData) 
    (v,T)=_valuevt(df)
    P= _pressure(df)
    return P*v/(Unitful.ustrip(Unitful.R)*T)
end


#
#Entropy
#

function _entropy(model::AbstractHelmholtzModel,v,T,x)
    return -_dfdt(model,v,T,x)
end

function _entropy(df::HelmholtzResultData)
    return -_gradientvt(df)[2]
end

function entropy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _entropy(model,v2,T2,x)*1.0u"J/mol/K"
end

function entropy(df::HelmholtzResultData)
    return _entropy(df)*1.0u"J/mol/K"
end

#
#enthalpy
#
function _enthalpy(model::AbstractHelmholtzModel,v,T,x)
    (f,df) = _gradientvt2(m,v,T,x)
    return f - df[1]*v-df[2]*T
end



function _enthalpy(df::HelmholtzResultData)
    x = _valuevt(df)
    a = _helmholtzvt(df)
    da= _gradientvt(df)
    return a - LinearAlgebra.dot(da,x)
end


function enthalpy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _enthalpy(model,v2,T2,x)*1.0u"J/mol"
end

function enthalpy(df::HelmholtzResultData)
    return _enthalpy(df::HelmholtzResultData)*1.0u"J/mol"
end

#
#Internal Energy
#
function _internal_energy(model::AbstractHelmholtzModel,v,T,x)
    u1 = core_helmholtz(model,v,ForwardDiff.Dual{typeof(T)}(T, one(T)),x)
    return core_helmholtz(model,v,T,x)-_gradientvt(model,v,T,x)[2]*T
end
#
#Implementation directly using dual numbers, is faster as expected, but dont derive it!
#
function _internal_energy2(model::AbstractHelmholtzModel,v,T,x)
    u1 = core_helmholtz(model,v,ForwardDiff.Dual{typeof(T)}(T, one(T)),x)
    return u1.value - T*u1.partials.values[1]
end

function _internal_energy(df::HelmholtzResultData)
    a = _helmholtzvt(df)
    (v,T)=_valuevt(df)
    da = _gradientvt(df)
    return a-v*da[2]
end

function internal_energy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _internal_energy2(model,v2,T2,x)*1.0u"J/mol"
end

function internal_energy(df::HelmholtzResultData) 
    return _internal_energy(df)*1.0u"J/mol"
end

##########
# Cv
##########
function _isochoric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    d2f = _hessianvt(model,v,T,x)
    return  -T*d2f[2,2]
    D
end

function _isochoric_heat_capacity(df::HelmholtzResultData)
    (v,T)= _valuevt(df)
    d2f = _hessianvt(df)
    return  -T*d2f[2,2]
    
end

function isochoric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return     _isocoric_heat_capacity(model,v2,T2,x)*1.0u"J/mol/K"
end

function isochoric_heat_capacity(df::HelmholtzResultData)
    return     _isocoric_heat_capacity(df)*1.0u"J/mol/K"
end
##########
# Cp,testing
##########
function _isobaric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    d2f = _hessianvt(model,v,T,x)
    return  -v*d2f[1,2]-T*d2f[2,2]  
end

function _isobaric_heat_capacity(df::HelmholtzResultData)
    (v,T)=_valuevt(df)
    d2f = _hessianvt(df)
    return  -v*d2f[1,2]-T*d2f[2,2]   
end

function isobaric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return     _isobaric_heat_capacity(model,v2,T2,x)*1.0u"J/mol/K"
end

function isobaric_heat_capacity(df::HelmholtzResultData)
  return _isobaric_heat_capacity(df::HelmholtzResultData)*1.0u"J/mol/K"
end
#######
# Speed of sound, testing
#######
function _sound_speed(model::AbstractHelmholtzModel,v,T,x)
    d2f = _hessianvt(model,v,T,x)
    return  v*sqrt(d2f[1,1]) 
end

function sound_speed(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return     _sound_speed(model,v2,T2,x)*1.0u"m/s"
end





####################################################

#mixing rules

####################################################

struct GeometricMeanRule <: AbstractMixingRule end
struct ArithmeticMeanRule <: AbstractMixingRule end
struct HarmonicMeanRule <: AbstractMixingRule end
struct PowerMeanRule{T} <: AbstractMixingRule where T<:Real
n::T
end #look at op_exponent below
#
#Abstract mixing rule, generates a mixing rule,based on 
# an operation, so mij = xi*xj*op(pi,pj)
#
function mixing_rule(op,x,p)
    N = length(x)
    @boundscheck checkbounds(p,N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1 : N
            res1 += p[i] * x[i]^2
                 for j = 1 : i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j])
            end
        end
    end
    return res1
end
#
#Abstract mixing rule, generates a mixing rule,based on 
# an operation, so mij = xi*xj*op(pi,pj)*Aij
#example: mixing_rule(geometric_mean_rule,x,Tc,1.-K)
#
function mixing_rule(op,x,p,A)
    N = length(x)
    checkbounds(A,N,N)
    @boundscheck checkbounds(p,N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1 : N
            res1 += p[i] * x[i]^2
                 for j = 1 : i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j])*A[i,j]
            end
        end
    end
    return res1
end
#
#Abstract asymetric mixing rule 
#Adds a Asymetric matrix A_asym, and a op_sim(xi,xj,Aij)
#the mayor example is the GERG2008 equation, where
#op_asym(xi,xj,Aij) = (xi+xj)/(Aij^2*xi + xj)
#
function mixing_rule_asymetric(op,op_asym,x,p,A,A_asym)
    N = length(x)
    checkbounds(A,N,N)
    checkbounds(A_asym,N,N)
    @boundscheck checkbounds(p,N)
    @inbounds begin
        res1 = zero(eltype(x))
        
        for i = 1 : N
            x[i] !=0 && begin
            res1 += p[i] * x[i]^2
                for j = 1 : i - 1
                    res1 += 2 * x[i] * x[j] * op(p[i], p[j])*A[i,j] *op_asym(x[i],x[j],A_asym[i,j])
                end
            end
        end
        end
    return res1
end


geometric_mean_rule(a,b)=sqrt(a*b)

arithmetic_mean_rule(a,b)=(a+b)/2

harmonic_mean_rule(a,b)=2*a*b/(a+b)

_power_mean_rule(a,b,n)=((a^(1/n) + b^(1/n))/2)^n


power_mean_rule(n) =(a,b)-> _power_mean_rule(a,b,n)

###############################
#Abstract mixing_matrix
#Creates A symetric matrix from an operation and a vector
#example: A[i,j]=sqrt(T[i]*T[j])
###############################
#
function mixing_matrix(op,p)
    N = length(p)
    A = Array{eltype(p),2}(undef,N,N)
    @boundscheck checkbounds(A,N,N)
    @inbounds begin
        res1 = zero(eltype(p))
        for i = 1 : N
            A[i,i] = p[i]
                for j = 1 : i - 1
                A[i,j] = op(p[i], p[j])
                A[j,i] = op(p[i], p[j])
            end
        end
    end
    return A
end

function mixing_matrix!(op,A,p)
    N = length(p)

    @boundscheck size(A)==(N,N)
    @inbounds begin
        res1 = zero(eltype(p))
        for i = 1 : N
            A[i,i] = p[i]
                for j = 1 : i - 1
                A[i,j] = op(p[i], p[j])
                A[j,i] = op(p[i], p[j])
            end
        end
    end
    return A

end






