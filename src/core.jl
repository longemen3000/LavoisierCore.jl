using JuliaDB, ForwardDiff, DiffResults, Unitful
include("utils.jl")

abstract type AbstractThermoModel end #base model
abstract type AbstractHelmholtzModel <:  AbstractThermoModel end
abstract type  AbstractGibbsModel <:  AbstractThermoModel end
abstract type  AbstractActivityModel <:  AbstractThermoModel end   

"""
core_helmholtz(model,V,T,x)

calculates the helmholtz energy. uses the following units:
v:molar volume, m3/mol ,SI
T:Temperature, K , SI
x:molar fraction, adimensional


"""
function core_helmholtz(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    throw(error("error x"))
    return -pi
end

function _gradient(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V T]
    f(d) = core_helmholtz(model,d[1],d[2],x)
    return ForwardDiff.gradient(f,d)
end

function _hessian(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V,T]
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    return ForwardDiff.hessian(f,d)
end

struct ResultDataHelmHoltz
    x::Array{Float64,1}
    fx::Float64
    dfdx::Array{Float64,1}
    d2fdx2::Array{Float64,2}
end

function diffdata(model::M,V,T,x::Array{R,1},order::Int64 = 2) where M<:AbstractHelmholtzModel where R <:Real
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
    return ResultDataHelmHoltz(d,f(d),DiffResults.gradient(res),DiffResults.hessian(res))
end

#
#Pressure
#
function _pressure(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return -_gradient(model,v,T,x)[1]
return 
end

function _pressure(diffdata::ResultDataHelmHoltz)
    return -diffdata.dfdx[1]
return 
end

function pressure(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R <: Real
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return uconvert(u"bar",_pressure(model,v2,T2,x)[1]*1.0u"Pa")
return 
end


#
#Entropy
#

function _entropy(model::M,V,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return -_gradient(model,V,T,x)[2]
end
function _entropy(diffdata::ResultDataHelmHoltz)
    return -diffdata.dfdx[2]
end
function entropy(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return _entropy(model,v2,T2,x)*1.0u"J/mol"
end



#
#enthalpy
#
function _enthalpy(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    df = _gradient(m,v,T,x)
    return core_helmholtz(m,v,T,x)-df[2]*T-v*df[1]
end

function _enthalpy(diffdata::ResultDataHelmHoltz)
    return diffdata.fx-diffdata.dfdx[2]*diffdata.x[2]-diffdata.dfdx[1]*diffdata.x[1]
end

function enthalpy(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return _enthalpy(model,v2,T2,x)*1.0u"J/mol"
end

#
#Internal Energy
#
function _internal_energy(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return core_helmholtz(m,v,T,x)-_gradient(m,v,T,x)[2]*T
end


function _internal_energy(diffdata::ResultDataHelmHoltz)
    return diffdata.fx-diffdata.x[2]*diffdata.dfdx[2] 
end

function internal_energy(model::M,v,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    (v2,T2) = _transformVT(v,T,[model.molecularWeight],x)
    return _internal_energy(model,v2,T2,x)*1.0u"J/mol"
end






