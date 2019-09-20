using JuliaDB, ForwardDiff, DiffResults, Unitful, LinearAlgebra, Zygote#master
include("utils.jl")

abstract type AbstractThermoModel end #base model
abstract type AbstractHelmholtzModel <:  AbstractThermoModel end
abstract type AbstractMixingRule end # for mixing rules
abstract type AbstractPhase end #for solvers
abstract type AbstractSpec{T,UNIT} end








"""
core_helmholtz(model,V,T,x)

calculates the helmholtz energy. uses the following units:
v:molar volume, m3/mol ,SI
T:Temperature, K , SI
x:molar fraction, adimensional


"""
function core_helmholtz(model::M,V,T,x) where M <:AbstractHelmholtzModel
end


Base.length(model::AbstractHelmholtzModel)=compounds_number(model)

#helmholtz but in total units (volume in m3, n in mol)
function  _helmholtzn(model::AbstractHelmholtzModel,V,T,n)
    sum_n = moles(mol_number(n))
    sum_n< zero(sum_n) && return zero(sum_n) 
    x = n*inv(sum_n)
    v = V*inv(sum_n) #this operation seems simple but it propagates the derivative of molar numbers,
    #since core_helmholtz is defined in the sense of molar fractions
    #with that in mind, the chemical potential of a mixture is n*gradient(_helmholtzn) in n
    return sum_n*core_helmholtz(model,v,T,x)
end

#helmholtz density, energy in a cubic meter, used in VTC equilibrium
function _helmholtzd(model::AbstractHelmholtzModel,T,C) 
    sum_c = sum(C)
    sum_c < zero(sum_c) && return zero(sum_c) 
    return sum_c*core_helmholtz(model,inv(sum_c),T,C/sum_c)
end

#this has the neat property: core_helmholtz(model,v,T,x)= -PV - sum(Ni*μi)
function _chemical_potential(model::AbstractHelmholtzModel,v,T,x)
f = z-> _helmholtzn(model,v,T,z)
return ForwardDiff.gradient(f,x)
end

function _dpdv(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(z->core_pressure(model,z,T,x),v)
end

function _d2pdv2(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(z->_dpdv(model,z,T,x),v)
end


function core_grad_vt(model::AbstractHelmholtzModel,v,T,x) 
    d = [v,T]
    f = z-> core_helmholtz(model,z[1],z[2],x)
    return ForwardDiff.gradient(f,d)
end

function core_grad_x(model::AbstractHelmholtzModel,v,T,x) 
    f = z-> core_helmholtz(model,v,T,z)
    return ForwardDiff.gradient(f,x)
end

# core_fg functions return value and gradient in one step, reusing calculations
function core_fg_x(model::AbstractHelmholtzModel,v,T,x) 
    df = DiffResults.GradientResult(x)
    df = ForwardDiff.gradient!(df,z->core_helmholtz(model,v,T,z),x)
    return (DiffResults.value(df),DiffResults.gradient(df))
end

function core_fg_vt(model::AbstractHelmholtzModel,v,T,x) 
    vt = [v,T]
    df = DiffResults.GradientResult(vt)
    df = ForwardDiff.gradient!(df,z->core_helmholtz(model,z[1],z[2],x),vt)
    return (DiffResults.value(df),DiffResults.gradient(df))
end

function core_dfdv(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.gradient(
        z->core_helmholtz(model,z[1],T,x),[v])[1]
end

function core_dfdt(model::AbstractHelmholtzModel,v,T,x)
    return ForwardDiff.derivative(
        z->core_helmholtz(model,v,z,x),T)
end

function core_hessian_vt(model::AbstractHelmholtzModel,v,T,x)
    f = z->core_helmholtz(model,z[1],z[2],x)
    return ForwardDiff.hessian(f,[v,T])
end

function core_fgh_vt(model::AbstractHelmholtzModel,v,T,x)
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
    f = z->  core_helmholtz(model,z[1],z[2],z[3:end])
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
(v2,T2)=_transformVT(v,T,model.molecular_weight,x)
return _diffdata(model,v2,T2,x,order)
end

_valuevt(df::HelmholtzResultData)  = df.value
_helmholtzvt(df::HelmholtzResultData)  = DiffResults.value(df.diffresult)
core_grad_vt(df::HelmholtzResultData)  = DiffResults.gradient(df.diffresult)[1:2]
core_hessian_vt(df::HelmholtzResultData)  = DiffResults.hessian(df.diffresult)[1:2,1:2]

#
#Pressure, a single derivative is faster
#
function core_pressure(model::AbstractHelmholtzModel,v,T,x)
    return -core_grad_vt(model,v,T,x)[1]
return 
end




function pressure(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return core_pressure(model,v2,T2,x)*1.0u"Pa"
end

function pressure(df::HelmholtzResultData) 
    return core_pressure(df)*1.0u"Pa"
end

function compressibility_factor(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return _pressure_impl(model,v2,T2,x)*v2/(ustrip(Unitful.R)*T2)
end

function compressibility_factor(df::HelmholtzResultData) 
    (v,T)=_valuevt(df)
    P= core_pressure(df)
    return P*v/(Unitful.ustrip(Unitful.R)*T)
end


#
#Entropy
#

function compressibility_factor(model::AbstractHelmholtzModel,v,T,x)
    return -core_grad_vt(model,v,T,x)[2]
end

function core_entropy(df::HelmholtzResultData)
    return -core_grad_vt(df)[2]
end

function entropy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return core_entropy(model,v2,T2,x)*1.0u"J/mol/K"
end

function entropy(df::HelmholtzResultData)
    return core_entropy(df)*1.0u"J/mol/K"
end

#
#enthalpy
#
function core_enthalpy(model::AbstractHelmholtzModel,v,T,x)
    (f,df) = core_fg_vt(model,v,T,x)
    return f - df[1]*v-df[2]*T
end



function core_enthalpy(df::HelmholtzResultData)
    x = _valuevt(df)
    a = _helmholtzvt(df)
    da= core_grad_vt(df)
    return a - LinearAlgebra.dot(da,x)
end


function enthalpy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return core_enthalpy(model,v2,T2,x)*1.0u"J/mol"
end

function enthalpy(df::HelmholtzResultData)
    return core_enthalpy(df::HelmholtzResultData)*1.0u"J/mol"
end

#
#Internal Energy
#
function core_internal_energy(model::AbstractHelmholtzModel,v,T,x)
    u1 = core_helmholtz(model,v,ForwardDiff.Dual{typeof(T)}(T, one(T)),x)
    return core_helmholtz(model,v,T,x)-core_grad_vt(model,v,T,x)[2]*T
end


function core_internal_energy(df::HelmholtzResultData)
    a = _helmholtzvt(df)
    (v,T)=_valuevt(df)
    da = core_grad_vt(df)
    return a-v*da[2]
end

function internal_energy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,model.molecularWeight,x)
    return core_internal_energy(model,v2,T2,x)*1.0u"J/mol"
end

function internal_energy(df::HelmholtzResultData) 
    return core_internal_energy(df)*1.0u"J/mol"
end

##########
# Cv
##########
function core_isochoric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    d2f = core_hessian_vt(model,v,T,x)
    return  -T*d2f[2,2]
    D
end

function core_isochoric_heat_capacity(df::HelmholtzResultData)
    (v,T)= _valuevt(df)
    d2f = core_hessian_vt(df)
    return  -T*d2f[2,2]
    
end

function isochoric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return core_isochoric_heat_capacity(model,v2,T2,x)*1.0u"J/mol/K"
end

function isochoric_heat_capacity(df::HelmholtzResultData)
    return core_isochoric_heat_capacity(df)*1.0u"J/mol/K"
end
##########
# Cp,testing
##########
function core_isobaric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    d2f = core_hessian_vt(model,v,T,x)
    return  -v*d2f[1,2]-T*d2f[2,2]  
end

function core_isobaric_heat_capacity(df::HelmholtzResultData)
    (v,T)=_valuevt(df)
    d2f = core_hessian_vt(df)
    return  -v*d2f[1,2]-T*d2f[2,2]   
end

function isobaric_heat_capacity(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return     core_isobaric_heat_capacity(model,v2,T2,x)*1.0u"J/mol/K"
end


function isobaric_heat_capacity(df::HelmholtzResultData)
  return core_isobaric_heat_capacity(df::HelmholtzResultData)*1.0u"J/mol/K"
end
#######
# Speed of sound, testing, needs checking
#######
function core_sound_speed(model::AbstractHelmholtzModel,v,T,x)
    d2f = core_hessian_vt(model,v,T,x)
    cp = -v*d2f[1,2]-T*d2f[2,2]
    cv = -T*d2f[2,2]
    
    m = 0.001*sum(molecular_weight(model).*x)
    return  v*sqrt(d2f[1,1]*(cp/cv)/m)
end

function sound_speed(model::AbstractHelmholtzModel,v,T,x)
    (v2,T2)=_transformVT(v,T,model.molecularWeight,x)
    return     core_sound_speed(model,v2,T2,x)*1.0u"m/s"
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

cubic_mean_rule(a,b)=((a^(1/3) + b^(1/3))/2)^3

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


#mass and molar types for dispatch
#the default is store the fractions and the amount.
#for molar types, the amount is in moles (gmol, as used in the metric system)
#for mass types, the amount is in kg

abstract type AbstractMaterialVector{T} <: AbstractVector{T} end

Base.size(A::AbstractMaterialVector) = size(A.values) #size of the vector
Base.length(A::AbstractMaterialVector)=length(A.values)
Base.getindex(A::AbstractMaterialVector, I::Int) = Base.getindex(A.values,I)
Base.getindex(A::AbstractMaterialVector, I::Vararg{Int})  = Base.getindex(A.values,I...)
Base.iterate(S::AbstractMaterialVector, state=1) = Base.iterate(S.values)
Base.IndexStyle(::Type{AbstractMaterialVector}) = IndexLinear()

struct MolFraction{T} <: AbstractMaterialVector{T}
    values::Vector{T}
    amount::T
end

struct MassFraction{T} <: AbstractMaterialVector{T}
    values::Vector{T}
    amount::T
end

struct MolNumber{T} <: AbstractMaterialVector{T}
    values::Vector{T}
    amount::T
end
Base.getindex(A::MolNumber, I::Int) =Base.getindex(A.amount .* A.values,I)
Base.getindex(A::MolNumber, I::Vararg{Int}) =Base.getindex(A.amount .* A.values,I...)



struct MassNumber{T} <: AbstractMaterialVector{T}
    values::Vector{T}
    amount::T
end
Base.getindex(A::MassNumber, I::Int) = Base.getindex(A.amount .* A.values,I)
Base.getindex(A::MassNumber, I::Vararg{Int}) = Base.getindex(A.amount .* A.values,I...)


moles(x::Union{MolFraction,MolNumber}) = x.amount
moles(x::Union{MolFraction,MolNumber},mw) = x.amount
moles(mw,x::Union{MolFraction,MolNumber}) = x.amount
moles(x::Union{MassFraction,MassNumber},mw) =sum(x.amount*1e3 .* x.values ./ mw)
moles(x::Union{MassFraction,MassNumber},mw::AbstractHelmholtzModel) = moles(x,molecular_weight(mw))
moles(mw,x::Union{MassFraction,MassNumber}) =moles(x,mw)

mass(x::Union{MassFraction,MassNumber}) = x.amount
mass(x::Union{MassFraction,MassNumber},mw) = x.amount
mass(x::Union{MassFraction,MassNumber},mw) = x.amount
mass(x::Union{MolFraction,MolNumber},mw) =sum(x.amount*1e-3 .* x.values .* mw)
mass(x::Union{MolFraction,MolNumber},mw::AbstractHelmholtzModel) = mass(x,molecular_weight(mw))
mass(mw,x::Union{MolFraction,MolNumber}) =mass(x,mw)

function mol_fraction(x::T1) where T1 <: AbstractVector
    T = eltype(x/sum(x))
     return MolFraction{T}(x/sum(x),one(T))
end

function mol_fraction(x::Union{MassFraction{T},MassNumber{T}},mw) where T
    sum1 = sum(x.amount*1e3 .* x.values ./ mw)
    MolFraction{T}(x.amount * 1e3 .* x.values ./sum1 ./ mw,one(T))
end

mol_fraction(x::MolFraction) = x
mol_fraction(model::AbstractHelmholtzModel,x) = mol_fraction(x)
mol_fraction(x,model::AbstractHelmholtzModel) = mol_fraction(x)
mol_fraction(x::MolNumber{T})  where T = MolFraction{T}(x.values,x.amount)
mol_fraction(x::Union{MolFraction,MolNumber},mw) = mol_fraction(x)
mol_fraction(mw,x::Union{MolFraction,MolNumber}) = mol_fraction(x)
mol_fraction(x::Union{MassFraction,MassNumber}) = throw(ArgumentError("molecular weight missing"))
mol_fraction(x::Union{MassFraction,MassNumber},model::AbstractHelmholtzModel) = mol_fraction(x,molecular_weight(model))
mol_fraction(mw,x::Union{MassFraction,MassNumber})=mol_fraction(x,mw)

function mol_number(x::T1)  where T1 <: AbstractVector
    T = eltype(x/sum(x))
     return MolNumber{T}(x/sum(x),sum(x))
end

function mol_number(x::Union{MassFraction{T},MassNumber{T}},mw) where T
    sum1 = sum(x.amount*1e3 .* x.values ./ mw)
    MolNumber{T}(x.amount * 1e3 .* x.values ./sum1 ./ mw,sum1)
end

mol_number(model::AbstractHelmholtzModel,x) = mol_number(x)
mol_number(x,model::AbstractHelmholtzModel) = mol_number(x)
mol_number(x::MolNumber) = x
mol_number(model::AbstractHelmholtzModel, x::Union{MolFraction, MolNumber}) = mol_number(x)
mol_number(x::MolFraction{T})  where T = MolNumber{T}(x.values,x.amount)
mol_number(x::Union{MolFraction,MolNumber},mw) = mol_number(x)
mol_number(mw,x::Union{MolFraction,MolNumber}) = mol_number(x)
mol_number(x::Union{MassFraction,MassNumber}) = throw(ArgumentError("molecular weight missing"))
mol_number(x::Union{MassFraction,MassNumber},model::AbstractHelmholtzModel) = mol_number(x,molecular_weight(model))
mol_number(mw,x::Union{MassFraction,MassNumber})=mol_number(x,mw)

function mass_fraction(x::T1) where T1 <: AbstractVector
    T = eltype(x/sum(x))
     return MassFraction{T}(x/sum(x),sum(x))
end

function mass_fraction(x::Union{MolFraction{T},MolNumber{T}},mw) where T
    sum1 = sum(x.amount*1e-3 .* x.values .* mw)
    MassFraction{T}(x.amount*1e-3 .* x.values ./sum1 .* mw,sum1)
end

mass_fraction(model::AbstractHelmholtzModel,x) = mass_fraction(x)
mass_fraction(x,model::AbstractHelmholtzModel) = mass_fraction(x)
mass_fraction(x::MassFraction) = x
mass_fraction(x::MassNumber{T})  where T = MassFraction{T}(x.values,x.amount)
mass_fraction(x::Union{MassFraction,MassNumber},mw) = mass_fraction(x)
mass_fraction(mw,x::Union{MassFraction,MassNumber}) = mass_fraction(x)
mass_fraction(x::Union{MolFraction,MolNumber}) = throw(ArgumentError("molecular weight missing"))
mass_fraction(x::Union{MolFraction,MolNumber},model::AbstractHelmholtzModel) = mass_fraction(x,molecular_weight(model))
mass_fraction(mw,x::Union{MolFraction,MolNumber})=mass_fraction(x,mw)

function mass_number(x::T1) where T1 <: AbstractVector
    T = eltype(x/sum(x))
     return MassNumber{T}(x/sum(x),sum(x))
end

function mass_number(x::Union{MolFraction{T},MolNumber{T}},mw) where T
    sum1 = sum(x.amount*1e-3 .* x.values .* mw)
    MassNumber{T}(x.amount*1e-3 .* x.values ./sum1 .* mw,sum1)
end

mass_number(model::AbstractHelmholtzModel,x) = mass_number(x)
mass_number(x,model::AbstractHelmholtzModel) = mass_number(x)
mass_number(x::MassNumber) = x
mass_number(x::MassFraction{T})  where T = MassFraction{T}(x.values,x.amount)
mass_number(x::Union{MassFraction,MassNumber},mw) = mass_number(x)
mass_number(mw,x::Union{MassFraction,MassNumber}) = mass_number(x)
mass_number(x::Union{MolFraction,MolNumber}) = throw(ArgumentError("molecular weight missing"))
mass_number(x::Union{MolFraction,MolNumber},model::AbstractHelmholtzModel) = mass_number(x,molecular_weight(model))
mass_number(mw,x::Union{MolFraction,MolNumber})=mass_number(x,mw)

####################################
#helmholtz phase
#is used to pass around results in a compact manner,maybe add helmholtz result data here?
#this completely defines a phase in helmholtz energy (v,T)
#temperature in Kelvin, molar volume in m3/mol,and the material is in mol number
struct HelmholtzPhase{T} <: AbstractPhase #this does nothing for now, what functions should implement an abstractPhase? 
    volume::T
    temperature::T
    material::MolNumber{T}
    molecularWeight::Array{T}
end


function helmholtz_phase(model,v,T,x)
TT = promote_type(eltype(x),typeof(v),typeof(T))
v1 = convert(TT,v)
T1 = convert(TT,T)
x1 = mol_number(TT.(x))
mw = TT.(molecular_weight(model))
return HelmholtzPhase{TT}(v1,T1,x1,mw)
end

function helmholtz_phase(model,v,T,x::AbstractMaterialVector{_T}) where _T
    TT = promote_type(_T,typeof(v),typeof(T))
    v1 = convert(TT,v)
    T1 = convert(TT,T)
    x1 = mol_number(model,x)
    mw = TT.(molecular_weight(model))
    return HelmholtzPhase{TT}(v1,T1,x1,mw)
    end

volume(a::HelmholtzPhase) = a.volume
temperature(a::HelmholtzPhase) = a.temperature
molecular_weight(a::HelmholtzPhase) = a.molecularWeight
moles(a::HelmholtzPhase)= moles(a.material)
mass(a::HelmholtzPhase) = mass(a.material,a.molecularWeight)

mol_fraction(a::HelmholtzPhase) = mol_fraction(a.material)
mol_number(a::HelmholtzPhase) = mol_number(a.material)
mass_fraction(a::HelmholtzPhase) = mass_fraction(a.material,molecular_weight(a))
mass_number(a::HelmholtzPhase) = mass_number(a.material,molecular_weight(a))
_unpack_helmholtz(a::HelmholtzPhase) = (volume(a),temperature(a),mol_fraction(a))

for op = (:_helmholtzn,    :core_helmholtz,    :core_grad_vt,    :core_grad_x,
    :core_fg_x,    :core_fg_vt,    :core_dfdv,    :core_dfdt,    :core_hessian_vt,
    :core_fgh_vt,    :core_pressure,    :compressibility_factor,    :core_entropy,
    :core_enthalpy,    :core_internal_energy,    :core_isochoric_heat_capacity,
    :core_isobaric_heat_capacity,    :core_sound_speed)
    @eval $op(model::AbstractHelmholtzModel,a::HelmholtzPhase) = $op(model,_unpack_helmholtz(a)...)
end
