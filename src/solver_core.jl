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
#volume solver 
#with P and T, finds a volume that satisfies pressure(v,T) = P
#without falling on the false loop.
#overload this for cubics and consistent EOS
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
function core_mol_volume(model::AbstractHelmholtzModel,P0,T0,x0)
    v = volume_solver(model,P0,T0,x0)
    return v
end

function core_mol_volume(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    v = volume_solver(model,P0,T0,x0,v0)
    return v
end
function core_mol_density(model::AbstractHelmholtzModel,P0,T0,x0)
    v = volume_solver(model,P0,T0,x0)
    return inv(v)
end

function core_mol_density(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    return inv(v)
end

function core_mass_volume(model::AbstractHelmholtzModel,P0,T0,x0)
    v = volume_solver(model,P0,T0,x0)
    k = dot(molecular_weight(model),x0)*1e-03 #kg/mol
    return v/k
end

function core_mass_volume(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    v = volume_solver(model,P0,T0,x0,v0)
    k = dot(molecular_weight(model),x0)*1e-03 #kg/mol
    return v/k
end

function core_mass_density(model::AbstractHelmholtzModel,P0,T0,x0)
    v = volume_solver(model,P0,T0,x0)
    k = dot(molecular_weight(model),x0)*1e-03 #kg/mol
    return k/v
end

function core_mass_density(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    v = volume_solver(model,P0,T0,x0,v0)
    k = dot(molecular_weight(model),x0)*1e-03
    return k/v
end


function mol_volume(model::AbstractHelmholtzModel,P0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    return core_mol_volume(model,P1,T1,x0)*1.0u"m^3/mol"
end

function mol_volume(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    v1 = _transform_v(v0,molecular_weight(model),x0)
    return core_mol_volume(model,P1,v1,T1,x0)*1.0u"m^3/mol"
end

function mol_density(model::AbstractHelmholtzModel,P0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    return core_mol_density(model,P1,T1,x0)*1.0u"mol/m^3"
end

function mol_density(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    v1 = _transform_v(v0,molecular_weight(model),x0)
    return core_mol_density(model,P1,v1,T1,x0)*1.0u"mol/m^3"
end

function mass_volume(model::AbstractHelmholtzModel,P0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    return core_mass_volume(model,P1,T1,x0)*1.0u"m^3/kg"
end

function mass_volume(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    v1 = _transform_v(v0,molecular_weight(model),x0)
    return core_mass_volume(model,P1,v1,T1,x0)*1.0u"m^3/kg"
end

function mass_density(model::AbstractHelmholtzModel,P0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    return core_mass_density(model,P1,T1,x0)*1.0u"kg/m^3"
end

function mass_density(model::AbstractHelmholtzModel,P0,v0,T0,x0)
    P1 = _transform_P(P0)
    T1 = _transform_T(T0)
    v1 = _transform_v(v0,molecular_weight(model),x0)
    return core_mass_density(model,P1,v1,T1,x0)*1.0u"kg/m^3"
end





function pt_flash(method::T,
    model::T2,
    P0,
    T0,
    x0,
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    throw(error("this method is not implemented."))
end

function pt_flash(
    model::T2,
    P0::Unitful.Pressure,
    T0::Unitful.Temperature,
    x0,
    method::T,
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    return core_pt_flash(method,
    model,
    _transform_P(P0),
    _transform_T(T0),
    x0,
    options)
end

function pt_flash(model::AbstractHelmholtzModel,
    phase0::HelmholtzPhase,
    method::AbstractHelmholtzSolver,options)
    P0 = core_pressure(model,phase0)
    T0 = temperature(phase0)
    x0 = mol_fraction(model,phase0)
    return core_pt_flash(model,P0,T0,x0,method,options)
end

#spinodal solver
#solves the equation dp/dv = 0 and gives the outer values as a tuple
function spinodal_solver(model::AbstractHelmholtzModel,T,x)
    fp(z) = core_pressure(model,z,T,x)
    dfp(z) = ForwardDiff.derivative(fp,z)
    min_v = dot(covolumes(model),x)
    max_v = 40*T0/100.0 # 0.01 atm

    #this is to be sure that (min_v,max_v) is a bracketing interval
    @show min_v
    @show max_v
    dv = find_zeros(fp,min_v,max_v,n_pts=21)
    dv
end


#function predict_saturation_pressure(model::AbstractHelmholtzModel,T0,x0,method)
function molar_volume_rackett(model::AbstractHelmholtzModel,T,x0)
    Pc = critical_pressure(model)
    Tc = critical_temperature(model)
    vc = critical_volume(model)
    Zc = core_compressibility_factor.(model,Pc,vc,Tc,Ref(1.0))
   # return vc
    return R_GAS .* Tc ./ Pc .* Zc .^ (1.0 .+ (1.0 .- T./ Tc).^(2/7) )
end