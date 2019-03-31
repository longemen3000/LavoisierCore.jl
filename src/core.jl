using JuliaDB, ForwardDiff, DiffResults
include("utils")

abstract type AbstractThermoModel end #base model
abstract type AbstractHelmholtzModel <:  AbstractThermoModel end
struct AbstractGibbsModel <:  AbstractThermoModel end
struct AbstractActivityModel <:  AbstractThermoModel end   

"""ischemicaldata(data)
Checks if a type can be used as a source to extract chemical data
    
""" 



function core_helmholtz(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    return 1+V+T
end

function _gradient(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V,T]
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    return ForwardDiff.gradient(f,d)
end

function _hessian(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    d = [V,T]
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    return ForwardDiff.hessian(f,d)
end

function _gradient(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    f(d) = core_helmholtz(model,d[1],d[2],d[3:end])
    return ForwardDiff.gradient(f,vcat(V,T,x))
end
function _der1(model::M,V,T,x::Array{R,1}) where M<:AbstractHelmholtzModel where R<:Real
    function f11(d)
        return core_helmholtz(model,d,T,x)
    end
    return ForwardDiff.derivative(f11,V)
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


function pressure(model::M,rho,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return -_gradient(model,rho,T,x)[1]*rho*rho
return 
end

function pressure(diffdata::ResultDataHelmHoltz)
    return -diffdata.dfdx[1]
return 
end

function entropy(model::M,V,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return -_gradient(model,V,T,x)[2]
end

function entropy(diffdata::ResultDataHelmHoltz)
    return -diffdata.dfdx[2]
return D
end



function internal_energy(model::M,V,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return core_helmholtz(model,V,T,x)-_gradient(model,V,T,x)[2]*T
end

function enthalpy(model::M,V,T,x::Array{R,1}) where M <: AbstractHelmholtzModel where R<:Real
    return core_helmholtz(model,V,T,x)-_gradient(model,V,T,x)[2]*T-_gradient(model,V,T,x)[1]*V
end


function internal_energy(diffdata::ResultDataHelmHoltz)
    return diffdata.fx-diffdata.x[2]*diffdata.dfdx[2]
return 
end




