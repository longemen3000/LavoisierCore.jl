# LavoisierCore.jl
Abstract Definitions and Functions to establish a thermodynamic library in julia

## What?
LavoisierCore intends to provide a thermodynamic property calculator, using helmholtz models and automatic differenciation to calculate all relevant properties with just one function definition. for example, lets say you are a chemist with a brand-new Helmholtz model,
but you don't want to implement all the properties that would transform your new EoS in something useful. 


## how?
for example, let's say you want to implement a new helmholtz model:
```
struct MyHelmholtzModel <: AbstractThermoModel
...#here you store your constants
...
end
```
then you

the package is the definition of the function `core_helmholtz(model,v,T,x)`, with v being the molar volume (m3/mol), T being the temperature (Kelvin), and x being the vector of molar fractions. the rest is obtained via autommatic differenciation.

```
core_helmholtz(model::MyHelmholtzModel,v,T,x)
...
end
```

And that's it!, with that, all relevant thermodynamic property's functions will be created, using the power of Julia's multiple dispatch and the powerful tools of automatic differenciation available (ForwardDiff at the moment, a reverse AD tool in the future).

at this moment, this package has implementations of 
* IAPWS 95 (formulation of water)
* GERG 2008 (21 natural gas compounds)

 and the following properties:

* compressibility factor
* pressure
* internal energy
* enthalpy
* entropy
* isochoric heat capacity
* isobaric heat capacity

if you defined your helmholtz equation, then you have access to all those functions.

## Usage example:
```
m = IAPWS95() #this model contains everything, so it doesn't need any variables to be created

# with unitful properties, it handles automatically the conversion beetween molar and mass density, 
#molar and mass especific volume, and the usual temperatures

julia> pressure(m,(322.0)u"kg/m^3",647.0u"K",[1.0])
220.64298608634772 bar

julia> pressure(m,55u"cm^3/mol",(647.0)u"K",[1.0])
220.59360203472647 bar

julia> pressure(m,55u"cm^3/mol",(373.85)u"°C",[1.0])
220.59360203472647 bar

#Big Floats Support
julia> pressure(m,big(322.0)u"kg/m^3",big(647.0)u"K",[1.0])
220.6429860863478813957299832315505132338714379480500058877684060654357069055415 bar

```
This package is in heavy development, don't dare to even try to use this in production, but feel free to submit your own Helmholtz equations,the long term idea is to implement the DiferencialEquations.jl of the thermodynamic equations

## Why?

Aspen is expensive, COOLPROP uses C++, and the thermodynamic papers are seldom implemented for free (for example,
the SAFT equations), and because is fun!, why not?

#News (September 19, 2019)

Added a equilibrium solver in Volume-Temperature-mol, based on a unified representation 
of helmholtz equilibria. a Pressure-Temperature-mol is following



