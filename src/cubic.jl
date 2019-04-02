include("core.jl")
include("utils.jl")
abstract type AbstractCubicModel <: AbstractHelmholtzModel end
abstract type AbstractAlphaT end
struct PengRobinson{Mixture,AlphaT} <: AbstractCubicModel
    criticalTemperature::Array{Float64,1}
    criticalPressure::Array{Float64,1}
    criticalVolume::Array{Float64,1}
    acentricFactor::Array{Float64,1}
    delta1::Float64
    delta2::Float64
    constanta::Float64
    constantb::Float64
    alphat::AlphaT
    mixture::Mixture
end


function PengRobinson(data,alphat::abstractAlphaT,mixture::AbstractMixture)
    Tc=disallowmissing(getproperty(data,:criticalTemperature))
    Pc=disallowmissing(getproperty(data,:criticalPressure))
    Vc=disallowmissing(getproperty(data,:criticalVolume))
    omega=disallowmissing(getproperty(data,:acentricFactor))
    return Cubic(Tc,Pc,Vc,omega,
    1+sqrt(2),1-sqrt(2),0.45724,0.07780,
    alphat,mixture)
end

struct Soave <: abstractAlphaT end
