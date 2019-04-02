include("core.jl")
include("IAPWS95.jl")
include("utils.jl")
using BenchmarkTools

m = IAPWS95()
rho0 = 838.025
T0 = 500
v0 = 1.0/molar_to_weight(rho0,[m.molecularWeight],[1.0])

Tx = 450

rhox = 890.341250 
vx = 1.0/molar_to_weight(rhox,[m.molecularWeight],[1.0])

hx = weight_to_molar(749161585.0,[m.molecularWeight],[1.0])

Px = 0.932203564 
rhot = 838.025
vt = weight_to_molar(1/rhot,[m.molecularWeight],[1.0])
rhot2 =1/molar_to_weight(vt,[m.molecularWeight],[1.0])

println(pressure(m,weight_to_molar(1/322,[m.molecularWeight],[1.0]),647,[1.0]))

hx = 1000*weight_to_molar(0.749161585e3  ,[m.molecularWeight],[1.0])
sx = 1000*weight_to_molar(0.210865845e1 ,[m.molecularWeight],[1.0])

testh() = @btime pressure(m,322.0u"kg/m^3",647.0u"K",[1.0])

