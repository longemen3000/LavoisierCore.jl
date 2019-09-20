using NLopt, ForwardDiff, DiffResults, Optim , Sobol, MappedArrays

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
_alpha_phase_2(α,x0,x1) = (x0-α*x1)/(1-α)
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

#molar volume and mol_number
function stvtd(model::AbstractHelmholtzModel,v,T,x0;verbose=false,test_phase=nothing,local_iters=nothing)      
   
    x = mol_fraction(x0,model)
    len = length(model)
    c0 =  x ./ v
    min_c = fill(1e-03,len)
    
    bi = covolumes(model)
    max_c = x./bi
    fn = c-> _helmholtzd(model,T,c)
    tpd = z -> _TPD(fn,c0,z)

    function _pressure(c) #pressure from molar concentrations
        v1 = inv(sum(c))
        x1 = c*v1
        return core_pressure(model,v1,T0,x1)
    end
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

    isnothing(local_iters) && (local_iters = 10*len)
   

    s = SobolSeq(min_c, max_c) #random distributed numbers with sobol sequences
    skip(s,local_iters) #recomended by the package
    crand = zeros(len)
    one_phase = true
    minf = Inf
    new_phase = zeros(len)
    
    for i in 1:local_iters
        next!(s,crand) #new random point
        minf,new_phase,ret = NLopt.optimize(opt_local, crand)       

        length(intersect(new_phase,min_c,max_c)) == 0 && (one_phase = false)
        minf <0.0                                     && (one_phase = false)
        (abs(minf)>sqrt(eps(Float64)))                && (one_phase = false)
        _pressure(new_phase) > 0.0                    && (one_phase = false) 
        
        if one_phase == false 
            if verbose == true
                println("Candidate phase found after $i local optimizations")
                println("Phase(moles/m3):: $new_phase")
                println("Potential Energy:: ",fn(new_phase))
            end
            break
        end
    end

    
    if one_phase == true
        if verbose == true
            println("Local optimizations failed to find a candidate phase,trying global optimization.")
        end
        opt_global = Opt(:GN_DIRECT_L, len) #sequencial quadratic programming
        opt_global.lower_bounds = min_c
        opt_global.upper_bounds = max_c
        opt_global.xtol_rel = 5e-5
        opt_global.min_objective = (x,grad) -> nlopt_form(tpd,x,grad,diffresults_cache_st)
        next!(s,crand)
        minf,new_phase,ret = NLopt.optimize(opt_local, crand)
        if length(intersect(new_phase,min_c,max_c)) == 0 &&
             minf <0.0 && 
             (abs(minf)>sqrt(eps(Float64))) &&
             _pressure(new_phase) > 0.0
            one_phase = false  
        end
    end
    
    if one_phase == false
        if verbose == true
            println("Candidate phase found in global optimization")
            println("Phase(moles/m3):: $new_phase")
        end
    else
        if verbose == true
            println("Global optimization failed to find a candidate phase,returning original phase.")
        end
        return helmholtz_phase(model,v,T,x) # the material is definitely one phase
    end

    βmin = [0.0]
    βmax = [min(1.0,minimum(c0./new_phase))]
    
    opt_eq0 = Opt(:GN_DIRECT_L, 1) #Global minimization in alpha
    opt_eq0.lower_bounds = βmin
    opt_eq0.upper_bounds = βmax
    opt_eq0.xtol_rel = 5e-2
    β0 = 0.5*(βmin + βmax)/2
    diffresults_cache_eq0 = DiffResults.GradientResult([1])
    fn_eq0 = β -> _st_twophases(fn,c0,new_phase,β[1])
    opt_eq0.min_objective = (x,grad) -> nlopt_form(fn_eq0,x,grad,diffresults_cache_eq0)
 
    minf_eq0,minβ,ret = NLopt.optimize(opt_eq0, β0)
    
    if !((minf_eq0<0.0) && 
        (!(minβ in [βmin,βmax,β0])) && 
        (ret != :FORCED_STOP) && 
        (abs(minf_eq0)>sqrt(eps(Float64))))
        
        if verbose == true
            println("Candidate phase does not meet the criteria to be an equilibrium phase")
            println("Phase(moles/m3):: $new_phase")
        end
        return helmholtz_phase(model,v,T,x)
    else
        if verbose == true
            println("Candidate phase meets the criteria to be an equilibrium phase")
            println("Proceeding to phase equilibria")
            println("Phase(moles/m3):: $new_phase")
            println("fraction:: ",minβ[1])
        end
    end

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
        fc2(x0,x1) = _alpha_phase_2(α,x0,x1)
        c2 = mappedarray(fc2, c0, c1)
        return α*fn(c1) + (1-α)*fn(c2)
    end

    function alpha_max(z)
        α = z[end];
        c1 = view(z,1:(length(z)-1))
        return α - minimum(c0./c1)
    end

    function covolume_max(z)
        #c1 is restricted via box bounds
        α = z[end];
        c1 = view(z,1:(length(z)-1))
        fc2(x0,x1) = _alpha_phase_2(α,x0,x1)
        c2 = mappedarray(fc2, c0, c1)
        c2[1]*bi[1]+c2[2]*bi[2] - 1 
        
    end

    
    min_eq = vcat(fill(0.0,len),sqrt(eps(Float64)))
    max_eq = vcat(fill(2.0^256,len),1.0-sqrt(eps(Float64)))
    diffresults_cache_eq = DiffResults.GradientResult(min_eq)
    diffresults_cache_const_eq = DiffResults.GradientResult(min_eq)

    opt_eq = Opt(:LD_SLSQP, len+1) #sequencial quadratic programming again, but with stricter tolerance
    opt_eq.lower_bounds = min_eq
    opt_eq.upper_bounds = max_eq
    opt_eq.ftol_abs = 1e-16
    objective_eq = (x,grad) -> nlopt_form(fn_eq,x,grad,diffresults_cache_eq)


    
    grad = zeros(len+1)

    exe = objective_eq(equilibrium_phase0,grad)
    println("c0:: ",equilibrium_phase0)
    println("function:: grad: ",grad,", f0: ",exe)
    println(minβ[1])
  
    opt_eq.min_objective = objective_eq
    inequality_constraint!(opt_eq,
        (x,grad) -> nlopt_form(alpha_max,x,grad,diffresults_cache_const_eq)
    )
    
    #phase equilibria
    minf_eq,equilibrium_phase,ret = NLopt.optimize(opt_eq,equilibrium_phase0)
    phase_success = true 

    function phase_to_helmholtz(c1,α)
        v1 = inv(sum(c1))
        n1 = mol_number(c1*v1*α)
        return helmholtz_phase(model,v1,T0,n1)
    end


    αfinal = equilibrium_phase[end]
    c1 = @view equilibrium_phase[1:end-1]
    fc2(x0,x1) = _alpha_phase_2(αfinal,x0,x1)
    c2 = mappedarray(fc2, c0, c1)
    
    v1 = inv(sum(c1))
    v2 = inv(sum(c2))
    isapprox(v1,v2;rtol = 1e-6) && (phase_success = false)
    ret != :FTOL_REACHED && (phase_success = false)
    #if the phase equilibria failed, return the original phase.
    phase_success = false && return helmholtz_phase(model,v,T,x)

    phase1  = phase_to_helmholtz(c1,αfinal)
    phase2  = phase_to_helmholtz(c2,1-αfinal)
    return [phase1,phase2]
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

function _pressure(res)
    x1 = res[1:2]
    v1 = inv(sum(res))
    x1 = x1*v1
    return pressure(minigerg,v1,T0,x1)
end





