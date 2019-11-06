abstract type AbstractThermoSolver end
abstract type AbstractHelmholtzSolver <: AbstractThermoSolver end
macro ifverbose(text)
println(typeof(text))
    if typeof(text) == String
    quote
        if verbose == true 
            println($text)
        end
    end |> esc
    else
        quote
        if verbose == true 
            @show $text
        end
        end |> esc
    end
end
function volume_solver(model::AbstractHelmholtzModel,P,T,x,v0 = nothing;verbose = false)  
    p(z) = core_pressure(model,z,T,x)
    fp(z) = p(z)-P
    dfp(z) = ForwardDiff.derivative(fp,z)
    if isnothing(v0)
    min_v = dot(covolumes(model),x)
    max_v = 40*T0/P0 #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while fp(max_v) > 0
        max_v *= 2
    end
    while fp(min_v) < 0
        min_v *=0.5
    end
    vv = find_zeros(fp,min_v,max_v)
    @ifverbose length(vv)
    if length(vv) == 1
        return vv[1]
    else length(vv) >= 1
        v1 = vv[1]
        v2 = vv[end]
        g1 = core_gibbs(model,P,v1,T,x)
        g2 = core_gibbs(model,P,v2,T,x)
        if g1 < g2
            return v2
        else
            return v1
        end
    end
else
    try 
        return Roots.newton(fp,dfp,v0)
    catch
        return volume_solver(model,P,T,x)  
    end
end
end

function pt_flash(method::T,
    model::T2,
    P0,
    T0,
    x0,
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    throw(error("this method is not implemented."))
end

function pt_flash(method::T,
    model::T2,
    P0::Unitful.Pressure,
    T0::Unitful.Temperature,
    x0,
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    return core_pt_flash(method,
    model,
    _transform_P(P0),
    _transform_T(T0),
    x0,
    options)
end

function pt_flash(method::AbstractHelmholtzSolver,
    model::AbstractHelmholtzModel,
    phase0::HelmholtzPhase,options)
    P0 = core_pressure(model,phase0)
    T0 = temperature(phase0)
    x0 = mol_fraction(model,phase0)
    return core_pt_flash(method,model,P0,T0,x0,options)
end

function vsol(model::AbstractHelmholtzModel,P,T,x,v0 = nothing;verbose = false)  
    p(z) = core_pressure(model,z,T,x)
    fp(z) = p(z)-P
    min_v = dot(covolumes(model),x)
    max_v = 40*T0/P0 #approx 5 times ideal gas
    

    return (fp,min_v,max_v)
end