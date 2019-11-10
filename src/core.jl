abstract type AbstractThermoModel end #base model
abstract type AbstractHelmholtzModel <:  AbstractThermoModel end
abstract type AbstractMixingRule end # for mixing rules
abstract type AbstractPhase end #for solvers

const R_GAS = ustrip(Unitful.R)
Base.broadcastable(x::AbstractHelmholtzModel) = Ref(x)




"""
core_helmholtz(model,V,T,x)

calculates the helmholtz energy. uses the following units:
v:molar volume, m3/mol ,SI
T:Temperature, K , SI
x:molar fraction, adimensional


"""


function core_helmholtz(model::M,V,T,x) where M <:AbstractHelmholtzModel
end
# partial_molar_property
#is equivalent to df/dni with with T,v, nj constant

function partial_molar_property(core_fn,model,v,T,x)
    function partial_fun(n)
        sum_n = sum(n)
        xx = n  .* inv(sum_n)
        return sum_n*core_fn(model,v/sum_n,T,xx) 
    end
    return ForwardDiff.gradient(partial_fun,x)
end

function core_logfugacity_coefficient(model,v,T,x)
    dAdn =  partial_molar_property(core_residual_helmholtz,model,v,T,x) 
    lnz = log(core_compressibility_factor(model,v,T,x))
    R = R_GAS
    return dAdn/(R*T) .- lnz
end


function core_logfugacity_coefficient(model,P,v,T,x)
    dAdn =  partial_molar_property(core_residual_helmholtz,model,v,T,x) 
    lnz = log(core_compressibility_factor(model,P,v,T,x))
    R = R_GAS
    return dAdn/(R*T) .- lnz
end

function core_logfugacity_coefficientn(model,v,T,n)
    return core_fugacity_coefficient(model,v,T,n/sum(n))
end

function core_logfugacity_coefficientn(model,P,v,T,n)
    return core_fugacity_coefficient(model,P,v,T,n/sum(n))
end

function core_logfugacity(model,v,T,x)
    dAdn =  partial_molar_property(core_residual_helmholtz,model,v,T,x) 
    R = R_GAS
    return log.(x .* (R*T/v)) .+ dAdn
end

function core_logfugacityn(model,v,T,n)
    return core_logfugacityn(model,v,T,n/sum(n))
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

function partial_molar_volume(model,v,T,x)
    return partial_molar_property(core_pressure,model,v,T,x)./core_dpdv(model,v,T,x)
end
#helmholtz density, energy in a cubic meter, used in VTC equilibrium
function _helmholtzd(model::AbstractHelmholtzModel,T,C) 
    sum_c = sum(C)
    #sum_c < zero(sum_c) && return zero(sum_c) 
    return sum_c*core_helmholtz(model,inv(sum_c),T,C/sum_c)
end

#gibbs density, energy in a cubic meter, used in VTC equilibrium
function _gibbsd(model::AbstractHelmholtzModel,T,C) 
    sum_c = sum(C)
    sum_c < zero(sum_c) && return zero(sum_c) 
    return sum_c*(core_helmholtz(model,inv(sum_c),T,C/sum_c)+core_pressure(model,inv(sum_c),T,C/sum_c))
end

function core_gibbs(model::AbstractHelmholtzModel,v,T,x)
    return core_helmholtz(model,v,T,x)+core_pressure(model,v,T,x)*v
end

function core_gibbs(model::AbstractHelmholtzModel,P,v,T,x)
    return core_helmholtz(model,v,T,x)+P*v
end
#this has the neat property: core_helmholtz(model,v,T,x)= -PV - sum(Ni*μi)
function _chemical_potential(model::AbstractHelmholtzModel,v,T,x)
f = z-> _helmholtzn(model,v,T,z)
return ForwardDiff.gradient(f,x)
end

function _chemical_potentiald(model::AbstractHelmholtzModel,v,T,x)
    f = z-> _helmholtzd(model,T,z/v)
    return ForwardDiff.gradient(f,x)
    end

function core_dpdv(model::AbstractHelmholtzModel,v,T,x)
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
    return ForwardDiff.derivative(z->core_helmholtz(model,z,T,x),v)
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
    return -core_dfdv(model,v,T,x)
return 
end




function pressure(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,molecular_weight(model),x)
    return core_pressure(model,v2,T2,x)*1.0u"Pa"
end

function pressure(df::HelmholtzResultData) 
    return core_pressure(df)*1.0u"Pa"
end

function core_compressibility_factor(model::AbstractHelmholtzModel,v,T,x) 
    return core_pressure(model,v,T,x)*v/(R_GAS*T)
end

function core_compressibility_factor(model::AbstractHelmholtzModel,P,v,T,x) 
    return P*v/(R_GAS*T)
end

function compressibility_factor(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,molecular_weight(model),x)
    return core_compressibility_factor(model,v2,T2,x)
end

function compressibility_factor(df::HelmholtzResultData) 
    (v,T)=_valuevt(df)
    P= core_pressure(df)
    return P*v/(Unitful.R_GAS*T)
end


#
#Entropy
#
function core_entropy(model::AbstractHelmholtzModel,v,T,x)
    return -core_dfdt(model,v,T,x)
return 
end
function core_entropy(df::HelmholtzResultData)
    return -core_grad_vt(df)[2]
end

function entropy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,molecular_weight(model),x)
    return core_entropy(model,v2,T2,x)*1.0u"J/mol/K"
end

function entropy(df::HelmholtzResultData)
    return core_entropy(df)*1.0u"J/mol/K"
end

#
#enthalpy
#
function core_enthalpy(model::AbstractHelmholtzModel,v,T,x)
    f = core_helmholtz(model,v,T,x)
    df = core_grad_vt(model,v,T,x)
    return f - df[1]*v-df[2]*T
end



function core_enthalpy(df::HelmholtzResultData)
    x = _valuevt(df)
    a = _helmholtzvt(df)
    da= core_grad_vt(df)
    return a - LinearAlgebra.dot(da,x)
end


function enthalpy(model::AbstractHelmholtzModel,v,T,x) 
    (v2,T2) = _transformVT(v,T,molecular_weight(model),x)
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
    (v2,T2) = _transformVT(v,T,molecular_weight(model),x)
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
    (v2,T2)=_transformVT(v,T,molecular_weight(model),x)
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
    (v2,T2)=_transformVT(v,T,molecular_weight(model),x)
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
    (v2,T2)=_transformVT(v,T,molecular_weight(model),x)
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

Base.size(A::AbstractMaterialVector) = size(A.n_values) #size of the vector
Base.length(A::AbstractMaterialVector)=length(A.n_values)
Base.IndexStyle(::Type{AbstractMaterialVector}) = IndexLinear()

struct MolFraction{T} <: AbstractMaterialVector{T}
    n_values::Vector{T}
    n_amount::T
    m_amount::T
    mw::Vector{T}
end

struct MassFraction{T} <: AbstractMaterialVector{T}
    n_values::Vector{T}
    n_amount::T
    m_amount::T
    mw::Vector{T}
end

struct MolNumber{T} <: AbstractMaterialVector{T}
    n_values::Vector{T}
    n_amount::T
    m_amount::T
    mw::Vector{T}
end

struct MassNumber{T} <: AbstractMaterialVector{T}
    n_values::Vector{T}
    n_amount::T
    m_amount::T
    mw::Vector{T}
end

Base.getindex(A::MolNumber, i::Int) =A.n_values[i]*A.n_amount
Base.getindex(A::MolNumber, I::Vararg{Int}) = A.n_values[I...]*A.n_amount
Base.getindex(A::MassNumber, i::Int) = A.n_values[i]*A.mw[i]*A.n_amount*1e-3
Base.getindex(A::MassNumber, I::Vararg{Int}) = A.n_values[I...]*A.mw[I...]*A.n_amount*1e-3
Base.getindex(A::MolFraction, i::Int) =A.n_values[i]
Base.getindex(A::MolFraction, I::Vararg{Int}) =A.n_values[I...]
Base.getindex(A::MassFraction, i::Int) = A.n_values[i]*A.mw[i]*A.n_amount*inv(A.m_amount)*1e-3
Base.getindex(A::MassFraction, I::Vararg{Int}) = A.n_values[I...]*A.mw[I...]*A.n_amount*inv(A.m_amount)*1e-3
core_moles(x::AbstractMaterialVector) = x.n_amount
core_mass(x::AbstractMaterialVector) = x.m_amount
moles(x::AbstractMaterialVector) = x.n_amount*1.0u"mol"
mass(x::AbstractMaterialVector) = x.m_amount*1.0u"kg"
#constructors
function mol_fraction(x::AbstractVector,mw::AbstractVector,n_amount_opt = nothing)
    T = promote_type(eltype(x),eltype(mw))
    n = convert(T,sum(x))
    n_values = convert.(T,x)
    n_values ./=n
    if isnothing(n_amount_opt)
        n_amount = sum(x)
    else
        n_amount = convert(T,n_amount_opt)
    end
    m_amount = LinearAlgebra.dot(x,mw)*1e-3
    return MolFraction{T}(n_values,n_amount,m_amount,convert.(T,mw))
end
function mol_number(x::AbstractVector,mw::AbstractVector,n_amount_opt = nothing)
    T = promote_type(eltype(x),eltype(mw))
    n = convert(T,sum(x))
    n_values = convert.(T,x)
    n_values ./=n
    if isnothing(n_amount_opt)
        n_amount = sum(x)
    else
        n_amount = convert(T,n_amount_opt)
    end
    m_amount = LinearAlgebra.dot(x,mw)*1e-3
    return MolNumber{T}(n_values,n_amount,m_amount,convert.(T,mw))
end

function mass_fraction(x::AbstractVector,mw::AbstractVector,m_amount_opt = nothing)
    T = promote_type(eltype(x),eltype(mw))
    if isnothing(m_amount_opt)
        m_amount = sum(x)
    else
        m_amount = convert(T,m_amount_opt)
    end
    n_values = normalizefrac(x .* inv.(mw))
    n_amount = LinearAlgebra.dot(x,inv.(mw))*1e-3
    return MassFraction{T}(n_values,n_amount,m_amount,convert.(T,mw))
end

function mass_number(x::AbstractVector,mw::AbstractVector,n_amount_opt = nothing)
    T = promote_type(eltype(x),eltype(mw))
    if isnothing(m_amount_opt)
        m_amount = sum(x)
    else
        m_amount = convert(T,m_amount_opt)
    end
    n_values = normalizefrac(x .* inv.(mw) .* 1e3)
    n_amount = LinearAlgebra.dot(x,inv.(mw))*1e3
    return MassNumber{T}(n_values,n_amount,m_amount,convert.(T,mw))
end

mass_number(x::AbstractVector,model::AbstractThermoModel) = mass_number(x,molecular_weight(model))
mass_fraction(x::AbstractVector,model::AbstractThermoModel) = mass_fraction(x,molecular_weight(model))
mol_number(x::AbstractVector,model::AbstractThermoModel) = mol_number(x,molecular_weight(model))
mol_fraction(x::AbstractVector,model::AbstractThermoModel) = mol_fraction(x,molecular_weight(model))

mass_number(model::AbstractThermoModel,x::AbstractVector) = mass_number(x,molecular_weight(model))
mass_fraction(model::AbstractThermoModel,x::AbstractVector) = mass_fraction(x,molecular_weight(model))
mol_number(model::AbstractThermoModel,x::AbstractVector) = mol_number(x,molecular_weight(model))
mol_fraction(model::AbstractThermoModel,x::AbstractVector) = mol_fraction(x,molecular_weight(model))

#conversions
function mass_number(x::AbstractMaterialVector{T}) where T
    return MassNumber{T}(x.n_values,x.n_amount,x.m_amount,x.mw)
end  
function mass_fraction(x::AbstractMaterialVector{T}) where T
    return MassFraction{T}(x.n_values,x.n_amount,x.m_amount,x.mw)
end
function mol_fraction(x::AbstractMaterialVector{T}) where T
    return MolFraction{T}(x.n_values,x.n_amount,x.m_amount,x.mw)
end
function mol_number(x::AbstractMaterialVector{T}) where T
    return MolNumber{T}(x.n_values,x.n_amount,x.m_amount,x.mw)
end

####################################
#helmholtz phase
#is used to pass around results in a compact manner,maybe add helmholtz result data here?
#this completely defines a phase in helmholtz energy (v,T)
#temperature in Kelvin, molar volume in m3/mol,and the material is in mol number
struct HelmholtzPhase{T} <: AbstractPhase #this does nothing for now, what functions should implement an abstractPhase? 
    volume::T
    temperature::T
    material :: AbstractMaterialVector{T}
end

function helmholtz_phase(model,v,T,x)
TT = promote_type(eltype(x),typeof(v),typeof(T))
v1 = convert(TT,v)
T1 = convert(TT,T)
x1 = mol_number(model,x)
return HelmholtzPhase{TT}(v1,T1,x1)
end

function helmholtz_phase(model,v,T,x::AbstractMaterialVector{_T}) where _T 
    TT = promote_type(_T,typeof(v),typeof(T))
    v1 = convert(TT,v)
    T1 = convert(TT,T)
    return HelmholtzPhase{TT}(v1,T1,mol_number(x))
end

core_moles(x) = sum(x) #if a normal vector comes to the solver, this can give you the moles
core_moles(a::HelmholtzPhase)= core_moles(a.material)
core_mass(a::HelmholtzPhase) = core_mass(a.material)
core_mol_volume(a::HelmholtzPhase) = a.volume
core_volume(a::HelmholtzPhase) = a.volume*core_moles(a)
core_mol_density(a::HelmholtzPhase) =inv(a.volume)
core_mass_volume(a::HelmholtzPhase) = a.volume*core_moles(a)/core_mass(a)
core_mass_density(a::HelmholtzPhase) =inv(a.volume*core_moles(a)/core_mass(a))
core_temperature(a::HelmholtzPhase) = a.temperature
core_molecular_weight(a::HelmholtzPhase) = a.material.mw

mol_volume(a::HelmholtzPhase) = core_molar_volume(a)*1.0u"mol/m^3"
mol_density(a::HelmholtzPhase) =core_molar_density(a)*1.0u"m^3/mol"
mass_volume(a::HelmholtzPhase) = core_mass_volume(a)*1.0u"m^3/kg"
mass_density(a::HelmholtzPhase) =core_mass_density(a)*1.0u"kg/m^3"
temperature(a::HelmholtzPhase) = core_temperature(a)*1.0u"K"
molecular_weight(a::HelmholtzPhase) = core_molecular_weight(a).*1.0u"g/mol"
moles(a::HelmholtzPhase)= core_moles(a)*1.0u"mol"
mass(a::HelmholtzPhase) = core_mass(a)*1.0u"kg"



mol_fraction(a::HelmholtzPhase) = mol_fraction(a.material)
mol_number(a::HelmholtzPhase) = mol_number(a.material)
mass_fraction(a::HelmholtzPhase) = mass_fraction(a.material)
mass_number(a::HelmholtzPhase) = mass_number(a.material)
_unpack_helmholtz(a::HelmholtzPhase) = (core_molar_volume(a),core_temperature(a),mol_fraction(a))

for op = (:_helmholtzn,    :core_helmholtz, :core_gibbs,    :core_grad_vt,    :core_grad_x,
    :core_fg_x,    :core_fg_vt,    :core_dfdv,    :core_dfdt,    :core_hessian_vt,
    :core_fgh_vt,    :core_pressure,    :compressibility_factor,    :core_entropy,
    :core_enthalpy,    :core_internal_energy,    :core_isochoric_heat_capacity,
    :core_isobaric_heat_capacity,    :core_sound_speed,
    :pressure,     :entropy,    :enthalpy,    :internal_energy,    :isochoric_heat_capacity,
    :isobaric_heat_capacity,    :sound_speed, :core_dpdv)
    @eval $op(model::AbstractHelmholtzModel,a::HelmholtzPhase) = $op(model,core_mol_volume(a),core_temperature(a),mol_fraction(a))
    @eval $op(model::AbstractHelmholtzModel,v,T,x::Union{MassFraction,MassFraction,MassNumber}) = $op(model,v,T,mol_fraction(x))
end