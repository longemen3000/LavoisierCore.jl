using NLopt,ForwardDiff,DiffResults

minigerg = GERG2008([19,1]) #H2S + Ch4
aire = GERG2008([2,16]) 
xaire = [0.79,0.21]
P0 = 4.53e6
T0 = 190
x0 = [0.05,0.95]
x01 = [1.731e-2,0.98269]
x02 = [0.06618,0.93382]

function fmin_2phases(x::Vector,model::AbstractHelmholtzModel,x0::Vector,T0)
    nt = length(x0)
    return _helmholtzn(model,x[nt+1],T0,x[1:(nt)])+_helmholtzn(model,x[end],T0,x[(nt+2):(end-1)]) -
    _helmholtzn(model,x0[end],T0,x0[1:end-1])
end

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





#uses x = [x V], #x0 is feed composition (one phase)
function sum_restriction(result::Vector,x::Vector,grad::Matrix,x0::Vector)
    nt = length(x0)
    Π = Int32(length(x)/nt)

    xzero = zero(eltype(x))
    xone = one(eltype(x))
    if length(grad) > 0
        grad .= ones(Π*nt,nt)
        for i = 1:nt
            sum_n = xzero
            for j = 1:Π
                nn = nt*(j-1)+i
                sum_n += x[nn]
            end
            result[i]=sum_n-x0[i]
        end
    else
        for i = 1:nt
            sum_n = xzero
            for j = 1:Π
                sum_n += x[nt*(j-1)+i]
            end
            result[i]=sum_n-x0[i]
        end
    end
    return result
end

function sum_restriction2(result,x,grad,x0)
    
    nt = length(x0)
    Π = Int32(length(x)/nt)

    xzero = zero(eltype(x))
    for i = 1:nt
        sum_n = xzero
        for j = 1:Π
            sum_n += x[nt*(j-1)+i]
        end
        result[i]=sum_n-x0[i]
    end
    
    if length(grad)>0
        grad=ones(size(grad))
    end
    return result
end

function _solver2(model,v0,T0,x0)
    r = rand()

    x00 = vcat(vcat(r*x0,r*v0),vcat((1-r)*x0,(1-r)*v0))
    nn = length(x00)
    lower = zeros(nn)
    upper = ones(nn)
    tol = upper*1e-06

    opt = NLopt.Opt(:LN_COBYLA, nn)
    opt.xtol_rel = 1e-6
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    equality_constraint!(opt, (r,x,g) -> sum_restriction2(r,x,g,x0),tol)
    fmin1 = (x,g) -> nlopt_form(vec -> fmin_2phases(vec,model,x0,T0),x,g)
    opt.min_objective = fmin1

    (min_f,min_x,status) = NLopt.optimize(opt, x00)

    return (min_f,min_x,status)
    return nn
end





