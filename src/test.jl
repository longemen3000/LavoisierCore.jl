using ForwardDiff 
using DiffResults 
#using Sobol
using MappedArrays
using LinearAlgebra
using Roots
using BenchmarkTools
using Unitful
#using Distributions
using PositiveFactorizations


include("utils.jl")
include("core.jl")
include("solver_core.jl")

include("IAPWS95.jl")
include("gerg2008.jl")
include("solver5.jl")

minigerg = GERG2008(:H2S,:CH4) #H2S + Ch4
P0 = 4.53e6
T0 = 190
x0 = [0.05,0.95]

x01 = [1.731e-02,0.98269]
x02 = [0.06618,0.93382]
v0 = 5.690353439153631e-5
v0 = 1.00367585015v0
v01 = 2.08e-4
v02 = 6.62e-05
#m = IAPWS95()
#gerg = GERG2008()
#rho0 = 838.025
#T0 = 500
#v0 = 1.0/molar_to_weight(rho0,[m.molecularWeight],[1.0])

#Tx = 450
#rhox = 890.341250 
#hx = 1000*weight_to_molar(0.749161585e3  ,[m.molecularWeight],[1.0])
#sx = 1000*weight_to_molar(0.210865845e1 ,[m.molecularWeight],[1.0])
#datx = (m,rhox*1.0u"kg/m^3",Tx*1.0u"K",[1.0])
#dat0 = (m,rho0*1.0u"kg/m^3",T0*1.0u"K",[1.0])



println("loaded")
