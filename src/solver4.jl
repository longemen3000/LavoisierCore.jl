using NLopt, ForwardDiff, DiffResults, Optim , Sobol, MappedArrays
include("ej1.jl")
#this completely defines a phase in helmholtz energy
#temperature in Kelvin, molar volume in m3/mol,and the material is in mol number
struct HelmholtzPhase{T}
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

volume(a::HelmholtzPhase) = a.volume
temperature(a::HelmholtzPhase) = a.temperature
molecular_weight(a::HelmholtzPhase) = a.molecularWeight
moles(a::HelmholtzPhase)= moles(a.material)
mass(a::HelmholtzPhase) = mass(a.material,a.molecularWeight)

mol_fraction(a::HelmholtzPhase) = mol_fraction(a.material)
mol_number(a::HelmholtzPhase) = mol_number(a.material)
mass_fraction(a::HelmholtzPhase) = mass_fraction(a.material,molecular_weight(a))
mass_number(a::HelmholtzPhase) = mass_number(a.material,molecular_weight(a))
_unpack_helmholtz(a::HelmholtzPhase) = (volume(a),temperature(a),molar_fraction(a))
#here i will procede to define an VTD solver:

#tangent plane distance functions of a potential
#temperature is constant, and the volume is reduced in the case of VTD to
#density = n/v, so the potencial is _helmholtzd(d)

#in TPD, 
_TPD(fn,x,y) = sum(ForwardDiff.gradient(fn,x) .* (x - y)) - fn(x) + fn(y)
_dTPDdy(fn,x,y) =   ForwardDiff.gradient(fn,y) .- ForwardDiff.gradient(fn,x)
_dTPDdy!(fn,F,dx,y) = begin F .=  ForwardDiff.gradient(fn,y) .- dx end
_st_twophases(fn,x,y,β) = β*fn(y) + (1-β)*fn((x .- β .* y)./ (1-β)) - fn(x)
_a_twophases(fn,x,y,β) = β*fn(y) + (1-β)*fn((x .- β .* y) ./ (1-β))

#auxiliary funtion: nlopt form to input functions with gradient
function nlopt_form(f,x,g,diffresult_cache)
    if length(g) > 0 
        ForwardDiff.gradient!(diffresult_cache,f,x)
        g .= DiffResults.gradient(diffresult_cache)
        return DiffResults.value(diffresult_cache)
    else
        return f(x)
    end
end

function autodiff_optim(f,x0)
    g! = (g,x) -> begin
        ForwardDiff.gradient!(g,f,x)
    end

    fg! = (g,x) -> begin
        ForwardDiff.gradient!(g,f,x)
        f(x)
    end
    h! = (H,x) -> begin
        ForwardDiff.hessian!(H,f,x)
    end 
    return TwiceDifferentiable(f, g!, fg!, h!,x0)
end

#molar volume and mol_fraction
function stvtd(model::AbstractHelmholtzModel,v,T,x0)      
    #println(core_helmholtz(model,v,T,x0))
    x = mol_fraction(x0,model)
    len = length(model)
    c0 =  x ./ v
    min_c = fill(1e-03,len)
    max_c = x ./ covolumes(model)
    fn = c-> _helmholtzd(model,T,c)
    tpd = z -> _TPD(fn,c0,z)

    #NLopt part. the strategy is the following:
    # a relatively low ammount of local optimizations is made, starting from random points.
    #if those optimizations fail to find a feasible point, a global optimization is performed.
    #if the global optimization fail to find a feasible point, there is not phase split.
    diffresults_cache_st = DiffResults.GradientResult(x)
    opt_local = Opt(:LD_SLSQP, len) #sequencial quadratic programming
    opt_local.lower_bounds = min_c
    opt_local.upper_bounds = max_c
    opt_local.xtol_rel = 1e-4
    opt_local.min_objective = (x,grad) -> nlopt_form(tpd,x,grad,diffresults_cache_st)

    local_iters = 5*len
    s = SobolSeq(min_c, max_c)
    skip(s,local_iters)
    crand = zeros(len)
    one_phase = true
    minf = Inf
    new_phase = zeros(len)
    for i in 1:local_iters
        #crand .= clamp.(c0 .* (2*rand(len)) .^ (1+i),min_c,max_c) #magic random sauce
        #trying Sobol
        next!(s,crand)
        minf,new_phase,ret = NLopt.optimize(opt_local, crand)
        
        if length(intersect(new_phase,min_c,max_c)) == 0 && minf <0.0 && (abs(minf)>sqrt(eps(Float64)))
            one_phase = false  
            break
        end
    end

    
    if one_phase == true
        
        opt_global = Opt(:GN_DIRECT_L, len) #sequencial quadratic programming
        opt_global.lower_bounds = min_c
        opt_global.upper_bounds = max_c
        opt_global.xtol_rel = 5e-5
        opt_global.min_objective = (x,grad) -> nlopt_form(tpd,x,grad,diffresults_cache_st)
        #crand .= clamp.(4 .* c0 .* rand(len) .^ 2,min_c,max_c) #magic random sauce
        next!(s,crand)
        minf,new_phase,ret = NLopt.optimize(opt_local, crand)
        if length(intersect(new_phase,min_c,max_c)) == 0 && minf <0.0 && (abs(minf)>sqrt(eps(Float64)))
            one_phase = false  
        end
    end
    
    
    if one_phase 
        return (_helmholtzn(model,v,T,x),helmholtz_phase(model,v,T,x)) # the material is definitely one phase
    end
    #println(crand,new_phase)
    βmin = [1e-06]
    βmax = [min(1.0,minimum(c0./new_phase))]
    
    
   diffresults_cache_eq0 = DiffResults.GradientResult([1])
    opt_eq0 = Opt(:GN_DIRECT_L, 1) #sequencial quadratic programming
    opt_eq0.lower_bounds = βmin
    opt_eq0.upper_bounds = βmax
    opt_eq0.xtol_rel = 5e-2
    β0 = 0.5*(βmin + βmax)/2
    fn_eq0 = β -> _st_twophases(fn,c0,new_phase,β[1])
    opt_eq0.min_objective = (x,grad) -> nlopt_form(fn_eq0,x,grad,diffresults_cache_eq0)
    minf_eq0,minβ,ret = NLopt.optimize(opt_eq0, β0)
    
    if !((minf_eq0<0.0) && 
        (!(minβ in [βmin,βmax,β0])) && 
        (ret != :FORCED_STOP) && 
        (abs(minf_eq0)>sqrt(eps(Float64))))
        return (_helmholtn(model,v,T,x),helmholtz_phase(model,v,T,x))
    end
    n2 = minβ[1] .* new_phase .* v
    n1 = x0 .- minβ[1] .* new_phase .* v
    #println(c0)
    #println(minβ[1]*new_phase)
    #println((mol_fraction(n2),mol_fraction(n1),minβ[1]))
    #println(_helmholtzn(model,minβ[1]v,T,minβ[1]*new_phase*v)+_helmholtzn(model,(1-minβ[1])*v,T,x0.-minβ[1]*new_phase*v))
    #return (minf,minβ[1],minβ[1]*new_phase)
    #here ends the stability test and the search of a candidade phase. now we proceed to calculate
    #the real phase split
    equilibrium_phase0 = [new_phase...,minβ[1]]
    
    
    
    function fn_eq(z)
        α = z[end];
        c1 = view(z,1:(length(z)-1))
        fc2(a,b) = (a-α*b)/(1-α)
        c2 = mappedarray(fc2, c0, c1)
        return α*fn(c1) + (1-α)*fn(c2)
    end
    function pressure_equality(z)
        α = z[end];
        c1 = view(z,1:(length(z)-1))
        fc2(a,b) = (a-α*b)/(1-α)
        c2 = mappedarray(fc2, c0, c1)
        v11 = inv(sum(c1))
        v22 = inv(sum(c2))
        return core_pressure(model,v11,T0,c1 .* v11) -  core_pressure(model,v22,T0,c2 .* v22)
    end
    println(min_c,minβ[1])
    min_eq = vcat(min_c,1e-6)
    max_eq = vcat(fill(2.0^256,len),1.0-1e-6)
    diffresults_cache_eq = DiffResults.GradientResult(min_eq)
    diffresults_cache_eq2 = DiffResults.GradientResult(min_eq)

    isoptim = false
    if isoptim == false
    opt_eq = Opt(:LD_SLSQP, len+1) #sequencial quadratic programming
    opt_eq.lower_bounds = min_eq
    opt_eq.upper_bounds = max_eq
    opt_eq.xtol_abs = 1e-16
    #println(fn_eq(equilibrium_phase0))
    opt_eq.min_objective = (x,grad) -> nlopt_form(fn_eq2,x,grad,diffresults_cache_eq)
    equality_constraint!(opt_eq, (x,grad) -> nlopt_form(pressure_equality,x,grad,diffresults_cache_eq2), 1e-8)
    println(fn_eq(equilibrium_phase0))
    println(pressure_equality(equilibrium_phase0))
    minf_eq,equilibrium_phase,ret = NLopt.optimize(opt_eq,equilibrium_phase0)
    #println(_helmholtzn(model,equilibrium_phase[end]*v,T,equilibrium_phase[end]*equilibrium_phase[1:end-1]*v)+_helmholtzn(model,(1-equilibrium_phase[end])*v,T,x0.-equilibrium_phase[end]*equilibrium_phase[1:end-1]*v))
    else
    #df = autodiff_optim(z->fn_eq(clamp.(z,min_eq,max_eq)),equilibrium_phase0)
    #dfc = TwiceDifferentiableConstraints(min_eq, max_eq)
    #res = Optim.optimize(df, equilibrium_phase0, NewtonTrustRegion())

    #bboptimize(fn_eq; SearchRange = ff = [box for box in zip(min_eq,max_eq)])

    prob = Objective(autodiff2(fn_eq),equilibrium_phase0)
    alg = MyNewton(zeros(len+1),0)
    H = zeros(len+1,len+1)
    g = zeros(len+1)
    
    return step!(prob,alg,H,g,x,false)

    end
    
    return (minf_eq,equilibrium_phase,ret)
end






function autodiff(f)
    res =(x,grad=[]) ->begin
        if length(grad) > 0 
            diffresult_cache = DiffResults.GradientResult(x) 
            ForwardDiff.gradient!(diffresult_cache,f,x)
            grad .= DiffResults.gradient(diffresult_cache)
            return DiffResults.value(diffresult_cache)
        else
            return f(x)
        end
    end
    return res 
    end





