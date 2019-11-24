abstract type AbstractThermoSolver end
abstract type AbstractHelmholtzSolver <: AbstractThermoSolver end
#volume solver 
#with P and T, finds a volume that satisfies pressure(v,T) = P
#without falling on the false loop.
#overload this for cubics and consistent EOS
function volume_solver(model::AbstractHelmholtzModel,P,T,x,v0 = nothing;no_pts = 7,all_volumes = false)  
    p(z) = core_pressure(model,z,T,x)
    fp(z) = p(z)-P
    dfp(z) = ForwardDiff.derivative(fp,z)
    if isnothing(v0)
    min_v = dot(covolumes(model),x)
    max_v = 40*T/P #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while fp(max_v) > 0
        max_v *= 2
    end
    while fp(min_v) < 0
        min_v *=0.5
    end
    vv = find_zeros(fp,min_v,max_v,no_pts=no_pts)
    if all_volumes == true 
        return vv
    end
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
        return find_zero(fp,v0)
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
function spinodal_solver(model::AbstractHelmholtzModel,T,x,predictor::Symbol=:PengRobinson)
    fp(z) = core_pressure(model,z,T,x)
    dfp(z) = core_dpdv(model,z,T,x)
    d2p(z) = ForwardDiff.derivative(dfp,z)
    min_v = dot(covolumes(model),x)
    P0 = dot(pressure_predictor(model,T,predictor),x) #by raoult, sum(xiP0i) = P
    max_v = 40*T/P0 # 0.01 atm
    while dfp(max_v) > 0
        max_v *= 2
    end
    while dfp(min_v) > 0
        min_v *=0.5
    end
    #this is to be sure that (min_v,max_v) is a bracketing interval
    dv = find_zeros(dfp,min_v,max_v,no_pts=21)
    if length(dv) > 2
        return (dv[1],dv[end],:success)
    elseif length(dv) == 0
        return (NaN,NaN,:zero_points)
    elseif length(dv) == 1
        return (dv[1],dv[1],:one_point)
    elseif length(dv) == 2
        if (d2p(dv[1]) > 0.0) && (d2p(dv[2]) < 0.0) #v1 checked
            return (dv[1],dv[2],:success)
        end
        return (dv[1],dv[2],:two_points)
    end
end

#pressure_predictor:
#tries to predict a saturation pressure given only the saturation temperature
#TODO: include a trait to let the model choose his own predictor?
#let the user build his own predictors?
#an ideal predictor would be an antoine equation, for example
pressure_predictor(model,T,pred::Symbol) = pressure_predictor(model,T,Val(pred))
function pressure_predictor(model,T,pred::Val{:PengRobinson})
    Tc = critical_temperature(model) 
    Tr = T ./ Tc
    m(ω) = 0.37464 + 1.54226*ω - 0.26992*ω^2  
    _α(ω,_Tr) = (1+m(ω)*(1-sqrt(_Tr)))^2
    α =_α.(acentric_factor(model),Tr)
    Pc = critical_pressure(model)
    aa =-2.605272488488440 
    bb = -9.017571450539830 
    cc = -19.896014683288000 
    dd = 0.579677284391001 
    ee = 22.501011849124900 
    ff = -0.041440165126830
    xx = α ./ Tr
    @. xx = (aa + cc*xx + ee*xx^2)/(1 + bb*xx + dd*xx^2 + ff*xx^3)
    xx .= exp.(xx) .* Tr
    return xx .* Pc
end

function core_pure_t_flash(model,T,x,options=nothing)
    if (Tc = dot(x,critical_temperature(model))) ≈ T
        return helmholtz_phase(model, dot(critical_volume(model),x), T, x)
    elseif T > Tc
        throw(error("the phase is supercritical at temperature = $T K"))
    end
    p_pred = dot(x,pressure_predictor(model,T,:PengRobinson))
    P = p_pred
    p(z) = core_pressure(model,z,T,x)
    fp(z) = p(z)-P
    min_v = dot(covolumes(model),x)
    max_v = 40*T/p_pred #approx 5 times ideal gas
    #this is to be sure that (min_v,max_v) is a bracketing interval
    while fp(max_v) > 0
        max_v *= 2
    end
    while fp(min_v) < 0
        min_v *=0.5
    end
    vv = find_zeros(fp,min_v,max_v,no_pts=21)
    #@show P
    if length(vv) <= 1 
        throw(error("the phase is stable at temperature = $T K"))
    elseif length(vv) == 2
        throw(error("proper initial points not found at temperature = $T K"))
    elseif length(vv) > 2
       v1 = vv[1]
       v2 = vv[end]
        p1 = p(v1)
        p2 = p(v2)
        px = p_pred
        A(z) = core_helmholtz(model,z,T,x)
        v1old = 0
        v2old = 0
        for i = 1:20  
            if i > 1
                v1old = v1
                v2old = v2
            end
            _v1 = Threads.@spawn volume_solver(model,px,T,x,0.9*v1)
            _v2 = Threads.@spawn volume_solver(model,px,T,x,1.1*v2)
            
            v1 = fetch(_v1)
            v2 = fetch(_v2)
            if abs(v1-v1old)/v1 < 1e-15 && i > 1
                #println("v1 condition")
                break
            elseif abs(v2-v2old)/v2 < 1e-15 && i > 1
                #println("v2 condition")
                break
            end
            pold = px
            px = (A(v1)-A(v2))/(v2-v1)
            #@show px
            if abs(pold-px) < 1e-10*px
                #println("P condition")
                break
            end

        end
        phase1 = helmholtz_phase(model,v1,T,x)
        phase2 = helmholtz_phase(model,v2,T,x)
        return [phase1,phase2]
    end
end

function core_pure_p_flash(model,P,x,options=nothing)
    if (Pc = dot(x,critical_pressure(model))) ≈ P
        return helmholtz_phase(model, dot(critical_volume(model),x), T, x)
    elseif P > Pc
        throw(error("the phase is supercritical at temperature = $T K"))
    end
    Tc = dot(x,critical_temperature(model))
    _t_pred0(T) =dot(x,pressure_predictor(model,T,:PengRobinson)) - P
    t_pred = find_zero(_t_pred0,Tc)
    vv = volume_solver(model,P,t_pred,x,all_volumes=true,no_pts=21)
    if length(vv) <= 1 
        throw(error("the phase is stable at temperature = $T K"))
    elseif length(vv) == 2
        throw(error("proper initial points not found at temperature = $T K"))
    elseif length(vv) > 2
       v1 = vv[1]
       v2 = vv[end]
        px = P
        _A(z,T) = core_helmholtz(model,z,T,x)
        v1old = 0.0
        v2old = 0.0
        Told = -100.0
        Tx = t_pred
        for i = 1:20 
            if i > 1
                v1old = v1
                v2old = v2
            end
            A = z-> _A(z,Tx)
            _v1 = Threads.@spawn volume_solver(model,P,Tx,x,0.9*v1)
            _v2 = Threads.@spawn volume_solver(model,P,Tx,x,1.1*v2)  
            v1 = fetch(_v1)
            v2 = fetch(_v2)
            if abs(v1-v1old)/v1 < 1e-15 && i > 1
                #println("v1 condition")
                break
            elseif abs(v2-v2old)/v2 < 1e-15 && i > 1
                #println("v2 condition")
                break
            end
            _px(T) = (_A(v1,T)-_A(v2,T))/(v2-v1) - P
            Told = Tx
            Tx = find_zero(_px,Tx)
            if abs(Told-Tx) < 1e-10*Tx
                #println("P condition")
                break
            end
        end
        
        phase1 = helmholtz_phase(model,v1,Tx,x)
        phase2 = helmholtz_phase(model,v2,Tx,x)
        return [phase1,phase2]
    end
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

function pure_t_flash(
    model::T2,
    T0::Unitful.Temperature,
    x0,
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    return core_pure_t_flash_2(
    model,
    _transform_T(T0),
    x0,
    options)
end

function pure_t_flash(
    model::T2,
    T0::Unitful.Temperature,
    x0 = [1.0],
    options=nothing) where T <: AbstractHelmholtzSolver where T2 <: AbstractHelmholtzModel
    return core_pure_t_flash_2(
    model,
    _transform_T(T0),
    x0,
    options)
end
