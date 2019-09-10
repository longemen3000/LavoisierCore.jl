#include("gerg2008.jl")
using NLopt,ForwardDiff,DiffResults, Plots, StatsBase, Roots
#held implementation of thermo solver, after that, a interfase to build custom solvers

abstract type AbstractStabilityTest end
abstract type AbstractThermoSolver end
abstract type AbstractHelmholtzSolver <: AbstractThermoSolver end


#xw = zeros(2)
#xw[1] = 1.0
#@btime pressure(m,(322.0)u"kg/m^3",647.0u"K",[1.0])

struct HELDSolver
eps_lambda::Float64
eps_composition::Float64
eps_b::Float64
eps_x::Float64
eps_v::Float64
eps_mu::Float64
eps_g::Float64
min_volume::Float64
max_volume::Float64
min_x::Array{Float64,1}
max_x::Array{Float64,1}
end

#volume from packing fraction
_vol(b,eta)=b/eta
_randvol(b) = _vol(b,rand())

#random molar fraction, 
function randfrac(N::Integer)
    x = rand(N)
    return x ./sum(x)
end

#random molar fraction, but considering the presence or absence of elements in a vector
function randfrac(x0::Array)
    length(x0)>0 && begin
    x = similar(x0)
    zerox = 0.0
    for i = 1:length(x0)
     x[i]=ifelse(x0[i]==zerox,zerox,rand())
    end
    x /= sum(x)
    return x 
end
end

function randpurefrac(x0::Array)
    x = zeros(length(x0))
    i = rand(LinearIndices(findall(x->x!=0.0,x0)))
    x[i] = 1.0
    return x
end

function randshiftfrac(x0::Array)
    x = x0 + 0.1*rand(length(x0))
    x /=sum(x)
    return x
end







function _stability_function(model::AbstractHelmholtzModel,P0,T0,x0,v,x)
    res = 0.0
    (A,dAdx)= _gradientx2(model,v,T0,x);
        res =A+P0*v
        for i = 1:length(dAdx)-1
            x0[i] != 0 && (res+=dAdx[i]*(x0[i]-x[i]))
        end
    return res
end


###here starts the optimization routine


_minvx(model::AbstractHelmholtzModel,P0,T0,x0) = vcat(min_volume(model,P0,T0,x0),(zeros(length(x0))))
_maxvx(model::AbstractHelmholtzModel,P0,T0,x0) =vcat(max_volume(model,P0,T0,x0),map(z->ifelse(z==0,0.0,1.0),x0))
_randvx(model::AbstractHelmholtzModel,P0,T0,x0) = vcat(random_volume(model,P0,T0,x0),randfrac(x0))

minigerg = GERG2008([19,1]) #H2S + Ch4
P0 = 4.53e6
T0 = 190
x0 = [0.05,0.95]
v0 = 5.690353439153631e-5
x01 = [1.731e-2,0.98269]
x02 = [0.06618,0.93382]

#gives a function with their grad, automatically, a la
function nlopt_form(f,xx,gg)
        if length(gg) > 0
            df = DiffResults.GradientResult(xx)
            df = ForwardDiff.gradient!(df,f,xx)
            gg .= DiffResults.gradient(df)
            return DiffResults.value(df)
        else
            return f(xx)
        end
end

function sumx(x::Vector, grad::Vector)
    if length(grad)>0
    grad[:] .= ones(size(grad))
    grad[1] = 0.0
    end
    return (sum(x[2:end])-1.0)
end

 



function st(model::AbstractHelmholtzModel,P0,T0,x0)
    
    fmin0 = vec -> _stability_function(model,P0,T0,x0,vec[1],vec[2:end])
    fmin1 = (x,g) -> nlopt_form(fmin0,x,g)
    
    stability_N = 50*(length(x0)+1)
    
    vx = _randvx(model,P0,T0,x0)
    lower = _minvx(model,P0,T0,x0)
    upper = _maxvx(model,P0,T0,x0)

    opt = NLopt.Opt(:LD_SLSQP, length(x0)+1)
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.xtol_rel = 1e-6
    opt.stopval = -1e-8
    opt.min_objective = fmin1
    NLopt.equality_constraint!(opt,sumx)
    best_f = Inf

    for i = 1:stability_N
        (min_f,min_vx,status) = NLopt.optimize(opt, vx)
        min_f < 0.0 && (return (min_f,min_vx,i))  
        vx[:] .= _randvx(model,P0,T0,x0)  
        #min_f<best_f && begin vx[1]= min_vx[1];best_f=min_f end
        #println((min_f,min_vx,status),i)   
    end

    return (best_f,vcat(0.0,x0),-1)
end

function st_test(model::AbstractHelmholtzModel,P0,T0,x0,n::Int=100)
    x = 0
    for i = 1:n
    a =  st(model,P0,T0,x0)[3]
    if a > 0
    x+=1
    end
end
    return x/n
    end


####################
##PHASE 2::searching candidate phases,
#variant with support for true zero values in the composition
#######################

function x_set(x0)
    nc = length(x0)
    nc2 = count(z->z>0,x0)
    Mmin = zeros(nc,nc2-1)
    Mmax = zeros(nc,nc2-1)
    set = findall(z->z>0,x0)
    for i = set[1:end-1]
        mmin_i = x0[i]/2.0
        mmax_i = mmin_i+0.5
        Mmin[set,i].=(1.0-mmin_i)/(nc2-1)
        Mmax[set,i].=(1.0-mmax_i)/(nc2-1)
        Mmin[i,i]=mmin_i
        Mmax[i,i]=mmax_i
        end
    return (hcat(Mmin,Mmax))
end

function _gibbs1v(model,P0,T0,x0)
    G = v -> core_helmholtz(model,v[1],T0,x0)+P0*v[1] 
    dG = (grad,v) -> begin 
        ForwardDiff.gradient!(grad,G,v)
    end
    lower = [min_volume(model,P0,T0,x0)]
    upper = [max_volume(model,P0,T0,x0)]
    
    opt = NLopt.Opt(:GN_CRS2_LM, 1)
    opt.xtol_rel = 5e-3
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    fmin1 = (x,g) -> nlopt_form(v -> core_helmholtz(model,v[1],T0,x0)+P0*v[1], x,g)
    opt.min_objective = fmin1

    opt2 = NLopt.Opt(:LD_SLSQP, 1)
    opt2.xtol_abs = 1e-8
    opt2.xtol_rel =  1e-8
    opt2.lower_bounds = lower
    opt2.upper_bounds = upper
    opt2.min_objective = fmin1
    

    initial_x = lower + 0.7*(upper-lower)
    (min_f,min_x) = NLopt.optimize(opt, initial_x)
    (min_f,min_x) = NLopt.optimize(opt2, min_x)

    return min_x[1]
end



###THIS NEEDS TO BE OPTIMIZED, change to NLopt soon
function m_set(model::AbstractHelmholtzModel,P0,T0,x0,vol=random_volume(model,P0,T0,x0))

    x_set0 = x_set(x0)
    v0 = zeros(size(x_set0,2))
    vmin0 = min_volume(model,P0,T0,x0)
    vmax0 = max_volume(model,P0,T0,x0)
    #v00 =_gibbs1v(model,P0,T0,x0)  
    for i = 1:size(x_set0,2)
        #vmin[i] = Roots.find_zero(z->_dgibbs_held1(model,z,P0,T0,Mmin[:,i]), (min_volume(model,P0,T0,x0), max_volume(model,P0,T0,x0)), Roots.A42()) #minimizing gibbs energy
        #vmax[i] = Roots.find_zero(z->_dgibbs_held1(model,z,P0,T0,Mmax[:,i]), (min_volume(model,P0,T0,x0), max_volume(model,P0,T0,x0)), Roots.A42())    
        #v0[i] = Roots.find_zero(z->_pressure(model,z,T0,x_set0[:,i])-P0, (vmin0,vmax0), Roots.A42())  
        v0[i] = vol     
        #vmin[i] = Roots.find_zero(z->_pressure(model,z,T0,Mmin[:,i])-P0, random_volume(model,P0,T0,x0))      #dont, work
        #vmax[i] = Roots.find_zero(z->_pressure(model,z,T0,Mmax[:,i])-P0, random_volume(model,P0,T0,x0))       #dont, work  
   
    end
    M = vcat(transpose(v0),x_set0)
    
    return M
end






##here using NLopt: is faster, but doest work, fix later



#base problem


function _gibbs1(model,P0,T0,x0)
    G = v -> core_helmholtz(model,v[1],T0,x0)+P0*v[1] 
    dG = (grad,v) -> begin 
        ForwardDiff.gradient!(grad,G,v)
    end
    lower = [min_volume(model,P0,T0,x0)]
    upper = [max_volume(model,P0,T0,x0)]
    
    opt = NLopt.Opt(:GN_CRS2_LM, 1)
    opt.xtol_rel = 1e-2
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    fmin1 = (x,g) -> nlopt_form(v -> core_helmholtz(model,v[1],T0,x0)+P0*v[1], x,g)
    opt.min_objective = fmin1

    opt2 = NLopt.Opt(:LD_SLSQP, 1)
    opt2.xtol_abs = 1e-8
    opt2.xtol_rel =  1e-8
    opt2.lower_bounds = lower
    opt2.upper_bounds = upper
    opt2.min_objective = fmin1
    

    initial_x = lower + 0.7*(upper-lower)
    (min_f,min_x) = NLopt.optimize(opt, initial_x)
    #(min_f,min_x) = NLopt.optimize(opt2, min_x)

    return (min_f,min_x[1])
end


## volume using gibbs minimization

#base problem
function L_v(model::AbstractHelmholtzModel,T0,P0,x0,v,x,lambda)
    res = core_helmholtz(model,v,T0,x) + P0*v
    for i = 1:length(x0)
    res+= lambda[i]*(x0[i]-x[i])
    end
    return res
end

function outer_problem(model::AbstractHelmholtzModel,P0,T0,x0,m_set,UBD,GP)
   
    opt = NLopt.Opt(:LD_SLSQP, length(x0)+1)
    fmin0 = vec -> vec[1]
    fmin1 = (x,g) -> nlopt_form(fmin0,x,g)
    opt.max_objective = fmin1
    bounds = ones(length(x0)+1)*1e16
    opt.lower_bounds = -bounds
    opt.upper_bounds = bounds
 
    constraint_gp = vec-> (vec[1]-GP)
    NLopt.inequality_constraint!(opt,(x,g) -> nlopt_form(constraint_gp,x,g),1e-8)
    
    for i in 1:size(m_set,2)
        constraint_i0 = vec-> (vec[1] - L_v(minigerg,T0,P0,x0,m_set[1,i],m_set[2:end,i],vec[2:end]))
        NLopt.inequality_constraint!(opt,(x,g) ->nlopt_form(constraint_i0,x,g),1e-8)
    end

    (min_f,min_vx,status) = NLopt.optimize(opt,vcat(UBD,-100*ones(length(x0))))
    return (min_f,min_vx,opt)
end


function inner_problem(model::AbstractHelmholtzModel,P0,T0,x0,UBD,lambda)
   
    opt = NLopt.Opt(:LD_SLSQP, length(x0)+1)
    fmin0 = vec -> L_v(model,P0,T0,x0,vec[1],vec[2:end],lambda)
    fmin1 = (x,g) -> nlopt_form(fmin0,x,g)
    opt.min_objective = fmin1
    
    opt.lower_bounds = _minvx(model,P0,T0,x0)
    opt.upper_bounds = _maxvx(model,P0,T0,x0)
 
    constraint_gp = vec -> L_v(model,P0,T0,x0,vec[1],vec[2:end],lambda) - UBD
    NLopt.inequality_constraint!(opt,(x,g) -> nlopt_form(constraint_gp,x,g),1e-8)
    NLopt.equality_constraint!(opt,sumx)
    #NLopt.equality_constraint!(opt,fmin1)

    (min_f,min_vx,status) = NLopt.optimize(opt,_randvx(model,P0,T0,x0))
    return (min_f,min_vx,status)
end
 
function held(model::AbstractHelmholtzModel,P0,T0,x0)

#step 1:: stability test    
step1 = st(minigerg,P0,T0,x0)
step1[3] == -1 && begin
(_,v)=_gibbs1(model,P0,T0,x0)
return (v,x0,1.0)
end
#step 2, inicialization
mayor_k = 0
M = m_set(model,P0,T0,x0)
(GP,_)=_gibbs1(model,P0,T0,x0)
#inicialize UBD, lambda
UBD = GP
LBD = 0.0
(vv,lambda,_)=outer_problem(model,P0,T0,x0,M,UBD,GP)
UBD = vv
(vv2,vx,_) = inner_problem(model,P0,T0,x0,UBD,lambda)
return(vx,lambda)
end
