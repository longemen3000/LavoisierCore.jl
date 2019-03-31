module LavoisierCore


include("core.jl")
include("IAPWS95.jl")
include("utils.jl")

export IAPWS95
export AbstractHelmholtzModel, pressure, core_helmholtz,diffdata,entropy

end






#end

#degreesoffredom(material::AbstractMaterial)
#material.state

#function eqsolve!(material::AbstractMaterial;kwargs)
#degreesofFreedom(material) == 0 && error("Non-zero degrees of freedom.")
#end

##function 
#function equilibrium(thermo::AbstractThermoModel;state::AbstractState;kwargs)
#end

##moistAir




#end # module

#idea:

# data = addcompounds(["CH4","O2","H2","CO2"])
#model = PengRobinson(data,alphaTSoave(),VdWMixtureRule())
#model = PengRobinson(data), if you dont want parameters, it shall have a sane default

#what about transport variables?, those are dependent of the temperature, pressure, etc,
#maybe using models
#but there are not dependent of the model
#X1 = [0.7,0.2,0.05,0.05]
#state1 = pvtxstate(T=400u"K",P=30u"bar",X1)

#material1 = newmaterial(state1,model)

#eqsolve!(material1)
#material2 = eqsolve(material1)
#entalphy(material1) #save cache of results in the form of the state vector and the value
#entropy(material2)  #probably using CAPE open definitions here


#another example, moist air, forking for psycro.jl to use his functions:

# data = addcompounds(["H2O","N2","O2"])
#MoistAirModel = ASHRAEMoistAir(data), #the system should check if the compounds are correct, it not, then throw an error
#this can be done in the database, for specific models
#state2 = moiststate(T=273u"K",wetBulb = 291.15u"K")
#moistair1 = newmaterial(model)

