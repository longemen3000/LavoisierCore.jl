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
return 1
end

function helmholtz(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return core_helmholtz(model,v2,T2,x)[1]*1.0u"J/mol"
end


#
#Stability test
#should return true if the phase is stable, false if not
#



function _gradient(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V T]
    f(d) = core_helmholtz(model,d[1],d[2],x)
    return ForwardDiff.gradient(f,d)
end

function _dfdv(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->core_helmholtz(model,z,T,x),v)
end

function _d2fdv2(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->_dfdv(model,z,T,x),v)
end

function _dfdt(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->core_helmholtz(model,v,z,x),T)
end

function _hessian(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V,T]
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    return ForwardDiff.hessian(f,d)
end


struct HelmholtzResultData
    value::Array{Float64,1}
    diffresult::DiffResults.DiffResult
end

function diffdata(model::AbstractHelmholtzModel ,V,T,x,order::Int64 = 2)
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

#
#Pressure
#
function _pressure(model::AbstractHelmholtzModel,v,T,x)
    return -_gradient(model,v,T,x)[1]
return 
end


function _pressure(diffdata::HelmholtzResultData)
    return -diffdata.dfdx[1]
return 
end

function pressure(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return uconvert(u"bar",-_dfdv(model,v2,T2,x)*1.0u"Pa")
return 
end




#
#Entropy
#

function _entropy(model::AbstractHelmholtzModel,v,T,x)
    return -_gradient(model,v,T,x)[2]
end

function _entropy(diffdata::HelmholtzResultData)
    return -diffdata.dfdx[2]
end

function entropy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return -_dfdt(model,v2,T2,x)*1.0u"J/mol"
end

#
#enthalpy
#
function _enthalpy(model::AbstractHelmholtzModel,v,T,x)
    df = _gradient(m,v,T,x)
    return core_helmholtz(m,v,T,x)-df[2]*T-v*df[1]
end

function _enthalpy(diffdata::HelmholtzResultData)
    return diffdata.fx-diffdata.dfdx[2]*diffdata.x[2]-diffdata.dfdx[1]*diffdata.x[1]
end


function enthalpy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return _enthalpy(model,v2,T2,x)*1.0u"J/mol"
end

#
#Internal Energy
#
function _internal_energy(model::AbstractHelmholtzModel,v,T,x)
    u1 = core_helmholtz(model,v,ForwardDiff.Dual{typeof(T)}(T, one(T)),x)
    return core_helmholtz(model,v,T,x)-_gradient(model,v,T,x)[2]*T
end
#
#Implementation directly using dual numbers, is faster as expected, but dont derive it!
#
function _internal_energy2(model::AbstractHelmholtzModel,v,T,x)
    u1 = core_helmholtz(model,v,ForwardDiff.Dual{typeof(T)}(T, one(T)),x)
    return u1.value - T*u1.partials.values[1]
end

function _internal_energy(diffdata::HelmholtzResultData)
    return diffdata.fx-diffdata.x[2]*diffdata.dfdx[2] 
end

function internal_energy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return _internal_energy2(model,v2,T2,x)*1.0u"J/mol"
end

#
#isobaric_heat_capacity
#

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
#
function mixing_rule_asymetric(op,op_asym,x,p,A,A_asym)
    N = length(x)
    checkbounds(A,N,N)
    checkbounds(A_asym,N,N)
    @boundscheck checkbounds(p,N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1 : N
            res1 += p[i] * x[i]^2
                 for j = 1 : i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j])*A[i,j] *op_asym(x[i],x[j],A_asym[i,j])
            end
        end
    end
    return res1
end


geometric_mean_rule(a,b)=sqrt(a*b)

arithmetic_mean_rule(a,b)=(a+b)/2

harmonic_mean_rule(a,b)=2*a*b/(a+b)

_power_mean_rule(a,b,n)=((a^n + b^n)/2)^(1/n)

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








