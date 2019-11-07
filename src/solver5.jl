struct Gernert <: AbstractHelmholtzSolver end

function core_pt_flash(method::Gernert,
    model::M,
    P0::Real,
    T0::Real,
    x0,
    options=nothing) where M <: AbstractHelmholtzModel
    return gernert_pt(model,P0,T0,x0)
end

function _rachford_rice(z,K) #find phase fraction, with z:feed phase and K:vaporization constant
    function f0(ψ)
        res = zero(eltype(K))
        for i = 1:length(K)
            if (K[i] == Inf) || (K[i] == -Inf)
                res += z[i]/ψ
            else
                res +=z[i]*(K[i]-1)/(1-ψ*(1-K[i]))
            end
        end
        return res
    end
    df0(ψ) = ForwardDiff.derivative(f0,ψ)
    try
        return Roots.find_zero((f0,df0),0.5,Roots.Newton()) 
    catch
        return Inf
    end
end  

function _rachford_rice(y,K,ψ)
    if  ψ == one(ψ)
        return sum(y .* (1 .- 1 ./ K))
    elseif ψ == zero(ψ)
        return sum(y .* (K .-1))
    else
        res = zero(eltype(K))
        for i = 1:length(K)
            if (K[i] == Inf) || (K[i] == -Inf)
                res += y[i]/ψ
            else
                res +=y[i]*(K[i]-1)/(1-ψ*(1-K[i]))
            end
        end
        return res
    end
end  



_K0(model,P0,T0,type::Symbol) = _K0(model,P0,T0,Val(type))
function _K0(model,P0,T0,type::Val{:mollerup})
    return (critical_pressure(model) ./ P0).*exp.(5.42.*(1.0 .- (critical_temperature(model)./T0)))
end
function _K0(model,P0,T0,type::Val{:wilson}) 
    return exp.( log.(critical_pressure(model) ./ P0)  .+ 5.373 .* (1.0 .+ acentric_factor(model)) .* (1.0 .- critical_temperature(model) ./ T0))
end

function gernert_pt(model::AbstractHelmholtzModel,P0::Real,T0::Real,x0)
    
    min_v = dot(x0,covolumes(model))
    max_v = 40*T0/P0
    v0 = Threads.@spawn volume_solver(model,P0,T0,x0)
    lnϕ0 = Threads.@spawn core_logfugacity_coefficient(model,fetch(v0),T0,x0)
    #=Step 1: Generation of trial phase compositions
    Starting with the assumption that two coexisting phases are present, 
    initial estimates of trial phase compositions are generated
    using the generalized Wilson correlation [25] to calculate K-values for all N components 
    =#
    K0 = _K0(model,P0,T0,:wilson) #this can be replaced if you have an starter that doesnt require critical data

    g_0 = _rachford_rice(x0,K0,0.0) 
    g_1 = _rachford_rice(x0,K0,1.0)  
    if g_0 <= 0  #bubble point assumption
        xi1 = copy(x0)
        xi2 = x0 .* K0 ./ (g_0 + 1)
    elseif g_1 >= 0 #dew point assumption
        xi1 = x0 ./ (-g_1 + 1) ./ K0   
        xi2 = copy(x0)
    else #two phase assumption
        ϕ = _rachford_rice(x0,K0)
        xi1 = x0 ./(1.0 .- ϕ .* (1.0 .- K0))
        xi2 = K0 .* x0 ./(1.0 .- ϕ .* (1.0 .- K0))
    end
    normalizefrac!(xi1)
    normalizefrac!(xi2)

    
    #= Step 2: Successive substitution method 
    Based on the initial estimates of the phase compositions, 
    three (3) steps of successive substitution are performed 
    in order to increase the accuracy of the estimates. 
    The successive substitution has been introduced 
    for the solution of phase equilibria conditions by Prausnitz and Chueh =#
    
    
    #initial step    
    vi1 = Threads.@spawn volume_solver(model,P0,T0,xi1)
    vi2 = Threads.@spawn volume_solver(model,P0,T0,xi2)

    lnϕi1 =similar(xi1)
    lnϕi2=similar(xi2)
    vv1 = fetch(vi1)
    vv2 = fetch(vi2)
    

    lnϕi1 = core_logfugacity_coefficient(model,P0,fetch(vv1),T0,xi1)
    lnϕi2 = core_logfugacity_coefficient(model,P0,fetch(vv2),T0,xi2)
    stable_phase = false
    
    for i = 1:3  
      
        K0 .= exp.(lnϕi1 .- lnϕi2)
     
        ϕ = _rachford_rice(x0,K0)
        if !(0 < ϕ < 1) #proceed to stability analysis
            stable_phase = true
            break
        end
        if i == 3 #partial loop
            break
        end
       
        xi1 .= x0 ./(1.0 .- ϕ .* (1.0 .- K0))
        xi2 .= K0 .* x0 ./(1.0 .- ϕ .* (1.0 .- K0))
        normalizefrac!(xi1)
        normalizefrac!(xi2)
        lnϕi1 .= core_logfugacity_coefficient(model,P0,fetch(vi1),T0,xi1)
        lnϕi2 .= core_logfugacity_coefficient(model,P0,fetch(vi2),T0,xi2)
    end
    v0 = fetch(v0)
    lnϕ0 = fetch(lnϕ0)
    @show lnϕ0
    function tpd(xi, lnϕ)
        return  sum(xi.*(log.(xi) .+lnϕ - log.(x0) .- lnϕ0))
    end
    if stable_phase == false
        tpd1 =   tpd(xi1,lnϕi1)
        tpd2 =   tpd(xi2,lnϕi2)
        ΔG = (1-ϕ)*tpd1 + ϕ*tpd2
        if (ΔG < 0) || (tpd1 < 0) || (tpd2 < 0)
            stable_phase = true
        end
    end
    #= the algorithm continues with a more detailed stability analysis. 
    In principle, the whole composition range needs to be checked 
    for negative tangent plane distances. 
    Since such a multi-dimensional search is impractical for multi-component mixtures 
    another more practical approach was suggested by Michelsen and Mollerup [23]. 
    From the initial Wilson K-values and the feed composition,
     heavy and light trial phase compositions are calculated =#
    
     vt_h = 0.0
     vt_l = 0.0
     vt_hmin = 0.0
     vt_lmin = 0.0
    if stable_phase == true
        K0=_K0(model,P0,T0,:wilson)
        xt_h = normalizefrac(x0./K0)
        xt_l = normalizefrac(x0 .*K0)
        xt_hmin = copy(xt_h)
        xt_lmin = copy(xt_l)
        vt_h = Threads.@spawn volume_solver(model,P0,T0,xt_h,vt_h)
        vt_l = Threads.@spawn volume_solver(model,P0,T0,xt_l,vt_l)
        
        vt_h_old = fetch(vt_h)
        vt_l_old = fetch(vt_l)
        vt_hmin =vt_h_old
        vt_lmin = vt_l_old
        lnϕt_l = core_logfugacity_coefficient(model,fetch(vt_l_old),T0,xt_l)
        lnϕt_h = core_logfugacity_coefficient(model,fetch(vt_h_old),T0,xt_h)
        tpd_hmin = tpd(xt_h,lnϕt_h)
        tpd_lmin = tpd(xt_l,lnϕt_l)
        xt_h .= x0.* exp.(lnϕ0) ./ exp.(lnϕt_h)
        xt_l .= x0.* exp.(lnϕ0) ./ exp.(lnϕt_l)
        normalizefrac!(xt_h)
        normalizefrac!(xt_l)
        i = 0
        iters = 5
        tpd_h = Inf
        tpd_l = Inf
        while i < iters
            i +=1     
            vt_h_old = fetch(vt_h) 
            vt_l_old = fetch(vt_l)        
            vt_h = Threads.@spawn volume_solver(model,P0,T0,xt_h,vt_h_old)
            vt_l = Threads.@spawn volume_solver(model,P0,T0,xt_l,vt_l_old)
            vt_h_old = fetch(vt_h) 
            vt_l_old = fetch(vt_l)  
            lnϕt_h = core_logfugacity_coefficient(model,fetch(vt_h_old),T0,xt_h)
            lnϕt_l = core_logfugacity_coefficient(model,fetch(vt_l_old),T0,xt_l)
            tpd_h = tpd(xt_h,lnϕt_h)
            tpd_l = tpd(xt_l,lnϕt_l)
            #println(tpd_h)
            #println(tpd_l)
            if tpd_h < tpd_hmin
                tpd_hmin = tpd_h
                xt_hmin .= xt_h
                vt_hmin = fetch(vt_h)
            end

            if tpd_l < tpd_lmin
                tpd_lmin = tpd_l
                xt_lmin .= xt_l
                vt_lmin = fetch(vt_l)
            end 
            if (tpd_lmin < 0) && (stable_phase == true)
                stable_phase = false
                iters +=3
            end
            if (tpd_hmin < 0) && (stable_phase == true)
                stable_phase = false
                iters +=3
            end 

            if xt_h ≈ x0 
                stable_phase = true
                break
            end
            if xt_l ≈ x0 
                stable_phase = true
                break
            end
            
            xt_h .= x0.* exp.(lnϕ0) ./ exp.(lnϕt_h)
            xt_l .= x0.* exp.(lnϕ0) ./ exp.(lnϕt_l)
           
            normalizefrac!(xt_h)
            normalizefrac!(xt_l)
            @show xt_h
            @show xt_l
            @show xt_lmin
            @show xt_hmin
        end
        
        vi1 = 0.0
        vi2 = 0.0
        if stable_phase == false
            xi1 .= xt_hmin
            xi2 .= xt_lmin
            vi1 = vt_hmin
            vi2 = vt_lmin
        end
    end

    if stable_phase == true
        return [helmholtz_phase(model,v0,T0,x0)] 
    end

    if stable_phase == false
       return _f0_isothermal6!(model,P0,T0,x0,xi1,xi2,vi1,vi2,fetch(v0))
    end
end

function _f0_isothermal6!(model,P0,T0,x0,x010,x020,v010,v020,v000)
    len = length(x0)
    x01 = zeros(eltype(x0),len)
    x02 = zeros(eltype(x0),len)
    ψ00 = Ref(zero(eltype(x0)))
    vv = zeros(eltype(x0),2)
    function f_0(z,v1,v2)
        x1 = vcat(z[1:end-1],1-sum(z[1:end-1]))
        β = z[end]
        @show x1
        return  β*core_gibbs(model,P0,v1,T0,x1) + (1-β)*core_gibbs(model,P0,v2,T0,(x0 .- β .* x1)/(1-β))
    end
    function fj!(F, J, z)
        β1 = z[end]
        sumx1 = one(eltype(x0))
        sumx2 = one(eltype(x0))
        for i = 1:len-1
            x01[i]= z[i]
            sumx1 -=x01[i]
            x02[i] =(x0[i] - β1*z[i])/(1- β1)
            sumx2 -=x02[i]
        end
        x01[end] = sumx1
        x02[end] = sumx2
        ψ00[]=β1
     
        vi1 = Threads.@spawn volume_solver(model,P0,T0,x01,vv[1])
        vi2 = Threads.@spawn volume_solver(model,P0,T0,x02,vv[2])
        v01 = fetch(vi1)
        v02 = fetch(vi2)
        vv[1] = v01
        vv[2] = v02
       
        #@show v01
        #@show v02
        
        f_01(z) = f_0(z,v01,v02)
        #@show f_01(z)
        if !(J == nothing)    
            ForwardDiff.hessian!(J,f_01,z)
        end
        if !(F == nothing)
            ForwardDiff.gradient!(F,f_01,z)    
        end
    end
F0 = similar(x0)
J0 = zeros(eltype(x0),len,len)
z0 = copy(x010)
ψ0 = abs((x0[1]-x020[1])/(x020[1]-x010[1]))
z0[end] = ψ0

dz = similar(z0)
cachedz = similar(z0)
fj!(F0,J0,z0)
#sol = nlsolve(only_fj(fj!),z0)
#
#sol)
τ = eps(eltype(x0))
ε2 = sqrt(τ)
ε3 = (τ)^(1/3)
ix = 0
isconverged = false
for i = 1:100
    fj!(F0,J0,z0)
    dz .= J0\-F0
    cholesky!(Positive, J0)
    α =1.0
    minz,maxz = extrema(z0 .+ α .* dz)
    while ! (0 <=minz < maxz <=1.0)
        α *=0.6321205588285577
        minz,maxz = extrema(z0 .+ α .* dz)
    end     
    @show z0
    @show vv
    f_j(z) = f_0(z,vv[1],vv[2])
    fj = f_j(z0)
    Θ = τ*(1+abs(fj)) 
    #@show Θ
    #@show  f_j(z0 .+ α .*dz)- fj
    #@show  fj
    #@show  f_j(z0 .+ α .*dz)
    #@show norm(dz)
    #@show sqrt(τ)*(1+norm(z0))
    
    if abs(f_j(z0 .+ α .*dz)- fj) < Θ
        isconverged =true
        break
    elseif norm(dz) < sqrt(τ)*(1+norm(z0))
        isconverged =true
        break
    elseif norm(x01-x02) < ε3
        #println("fractions are the same")
        break
    elseif abs(vv[1]-vv[2]) < ε2
        #println("volumes are the same")
        break 
    elseif i == 100
        break
    end
    ix +=1

    z0 .= z0 .+ α .*dz
end
    if isconverged == true
        n1 = mol_number(x01,molecular_weight(model),ψ00[])
        n2 = mol_number(x02,molecular_weight(model),1-ψ00[])
        return [helmholtz_phase(model,vv[1],T0,n1),helmholtz_phase(model,vv[2],T0,n2)] 
    else
        return [helmholtz_phase(model,v000,T0,x0)] 
    end
end

#mono_t_psiflash(model::AbstractHelmholtzModel,T0::Real,ψ0::Real,x0)
#Tr = T0/dot(x0,critical_temperature(model))
#if Tr>one(Tr)
#
#else
#end