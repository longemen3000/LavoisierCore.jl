
include("core.jl")
include("IAPWS95.jl")
include("gerg2008.jl")
include("held.jl")

using BenchmarkTools

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
