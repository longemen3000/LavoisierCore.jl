
include("core.jl")
include("IAPWS95.jl")
include("gerg2008.jl")
include("utils.jl")
using BenchmarkTools, Optim

#m = IAPWS95()
gerg = GERG2008()
rho0 = 838.025
T0 = 500
v0 = 1.0/molar_to_weight(rho0,[m.molecularWeight],[1.0])

Tx = 450
rhox = 890.341250 
hx = 1000*weight_to_molar(0.749161585e3  ,[m.molecularWeight],[1.0])
sx = 1000*weight_to_molar(0.210865845e1 ,[m.molecularWeight],[1.0])
datx = (m,rhox*1.0u"kg/m^3",Tx*1.0u"K",[1.0])
dat0 = (m,rho0*1.0u"kg/m^3",T0*1.0u"K",[1.0])

function _st(model::AbstractHelmholtzModel,P0,T0,x,vmin,vmax;pole=1e-3,iter = 10*length(x))
    lower = vcat(vmin,zeros(length(x)))
    upper = vcat(vmax,ones(length(x)))
    
    exp(pole/LinearAlgebra.norm(z))
    return true
end


xw = zeros(21)
xw[18]=1.0
#@btime pressure(m,(322.0)u"kg/m^3",647.0u"K",[1.0])
