# LavoisierCore.jl
Abstract Definitions and Functions to establish a thermodynamic library in julia

## What?
LavoisierCore intends to provide a thermodynamic property calculator, using helmholtz models and automatic differenciation to calculate
all relevant properties with just one function definition. for example, lets say you are a chemist with a brand-new Helmholtz model,
but you don't want to implement all the properties that would transform your new EoS in something useful. for example:
```
struct MyShinySAFT <: AbstractThermoModel
...
...
end
```

now, to use the power of this package, you just need to define one thing: your equation, of course. The principal interface to 
the package is the definition of the function `core_helmholtz(model,v,T,x)`, with v being the molar volume (m3/mol at the moment)
T being the temperature (Kelvin), and x being the vector of molar fractions.

```
core_helmholtz(model::MyShinySAFT,v,T,x)
...
end
```

And that's it!, with that, all relevant thermodynamic property's functions will be created, using the power of Julia's multiple dispatch
and the powerful tools of automatic differenciation available (ForwardDiff at the moment, a reverse AD tool in the future).

This package is in heavy development, don't dare to even try to use this in production, but feel free to submit your own Helmholtz
equations, a dream i have is to implement the DiferencialEquations.jl of the thermodynamic equations

## Why?

Aspen is expensive, COOLPROP uses C++, and the thermodynamic papers are seldom implemented for free (for example,
the SAFT equations), and because is fun!, why not?





