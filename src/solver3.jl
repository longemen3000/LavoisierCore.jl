using Optim, ForwardDiff, DiffResults, NLopt, NLsolve, LineSearches, PositiveFactorizations
using LinearAlgebra, BenchmarkTools, DiffResults, Roots
abstract type AbstractEquilibriumSpec end

struct AbstractHelmholtzPhases{N,V,T,W}
    phases::Int64
    molar_fractions::Vector{Vector{N}}
    volumes::Vector{V}
    temperature::T
    partition::Vector{W}
    #function AbstractHelmholtzPhases(v,T,x,phi)
     #   phas  = length(v)
     #   return new(phas,v,T,x,phi)
    #nd
end

#nn is the number of phases minus one
function _build_correction_matrix(nn)
     
    AA =[Diagonal(ones(nn)) Diagonal(zeros(nn))]
    BB =[Transpose(fill(10.0,nn)) Transpose(zeros(nn))]
    CC = [Diagonal(zeros(nn)) Diagonal(ones(nn))]
    DD = hcat(Transpose(zeros(nn)),Transpose(-ones(nn)))
    return dropzeros!(vcat(AA,BB,CC,DD))
end

#basis for stability testing
#TPD2(fn,x,y) = sum(ForwardDiff.gradient(fn,x),(x - y)) - fn(x) - fn(y)
TPD_vt(fn,x,y) = sum(ForwardDiff.gradient(fn,x) .* (x - y)) - fn(x) + fn(y)

dTPDdy_vt(fn,x,y) =   ForwardDiff.gradient(fn,y) .- ForwardDiff.gradient(fn,x)

function dTPDdy_vt!(fn,F,dx,y)
    F .=  ForwardDiff.gradient(fn,y) - dx
end
d2TPDdy2_vt(fn,x,y) = ForwardDiff.hessian(fn,y)

function d2TPDdy2_vt!(fn,J,x,y)
    J .= ForwardDiff.hessian(fn,y)
    #cholesky!(Positive,J)
end

function TPD_allderivatives_vt!(fn,diffresult,F,J,x,y,dx)
    ForwardDiff.hessian!(diffresult,fn,y)
    F .= DiffResults.gradient(diffresult) - dx
    J .= DiffResults.hessian(diffresult)
end





fn_vt(C,model::AbstractHelmholtzModel,T) = begin 
    C2= max.(0.0,C)
    _helmholtzd(model,T,C2)   
end
    

#stability in two phases
function st_vt(model::AbstractHelmholtzModel,V0,T0,n0)
    b0 = covolumes(model)
    n = length(n0) #length of compounds
    q = 0 #counter of random aproximations
    s = 10*n #maximum random aproximations
    k_max = 10 #max inner iterations of cholesky minimizations
    ε = 1e-06
    diff_c = 0
    yy = zeros(n) #holder of solutions
    delta = zeros(n) #holder of delta
    xx0 = zeros(n) #holder of x0
    randx0 = zeros(n) #holder of random vector
    F0 = zeros(n) #holder of function gradient
    J0 = zeros(n,n) #holder of function hessian
    yyplus1 = zeros(n) #another holder
    diffresult= DiffResults.HessianResult(yyplus1)
    
    for i in 1:n
        xx0[i] = n0[i]/V0
    end
    
    #xx0[end] = 1.0
    fn = z-> fn_vt(z,model,T0)
    dx = ForwardDiff.gradient(fn,xx0)
    tpd = z ->TPD_vt(fn,xx0,z)
    tpd0 = sum(ForwardDiff.gradient(fn,xx0) .* xx0) - fn(xx0)
    f! = (F,z) -> dTPDdy_vt!(fn,F,dx,z)
    j! = (J,z) -> d2TPDdy2_vt!(fn,J,xx0,z)
    fj! = (F,J,z) -> TPD_allderivatives_vt!(fn,diffresult,F,J,xx0,z,dx)
    norm_sol = z -> norm_vt(xx0,z)
    y_c = 0
    while true
        #building q-random aproximation to n0
        q==s && break
        q+=1
        k = 0
        
        randx0 .=4*rand(n).^2
        yy .= randx0 .* xx0
        

        normq = 1.0
        iterq = 0


        while true
            k == k_max && break
            k+=1

            fj!(F0,J0,yy)
            diff_c +=1
            cholesky!(Positive,J0)
            delta .= J0\(-F0)
            #println(yy)
            #println(delta)
            #println(yy+delta)
            
            #initial λ0, so maximum(yy+λ0*delta)>0
            λ0 = 1.0
            for i = 1:n
                if yy[i]+λ0*delta[i] <0
                    λ0 = (-yy[i]/delta[i])
                end
            end
            #println(λ0)
            #println("")
            λ=λ0
            #println("lambda0: ",λ0)
            iter_λ = 0
            tpdk0 = tpd(yy)
            tpdk1 = tpd(yyplus1)
            yyplus1 .= yy+λ*delta
            
            while iter_λ < 20
                iter_λ +=1                
                λ*=0.5
                tpdk1 < tpdk0 && break
                if tpdk1 < 0 
                    if (all(yyplus1.>sqrt(eps(eltype(n0)))/V0) && sum(b0.*yyplus1) <1)
                      
                     return (yyplus1,2,q)
                    end
                end
                tpdk1 = tpd(yyplus1)
                yyplus1 .= yy+λ*delta
            end

            tol = sqrt(sum(λ*delta./yy).^2)
            tol< 1e-06 && break
            yy .=yyplus1 #new value of yy
    
        end #end of qth aproximation
        
        if tpd(yy) < 0.
            if (all(yy.>sqrt(eps(eltype(n0)))/V0) && sum(b0.*yy) <1)
                return (yy,2,q)
            end
            
        end       
    end
    return(xx0,1,s,-1)
end #end of all aproximations




function test_st(model::AbstractHelmholtzModel,v,T,n)
    x = n/v
    res = st_vt(model,v,T,n)
    res[2] == 1 && return 0
    y  = res[1]
    βmax = min(1.0,minimum(x./y))
    
    f0β = β -> β*_helmholtzd(model,T,max.(0.0,y)) +
    (1-β)*_helmholtzd(model,T,max.(0.0,((x .- β*y)/((1-β))))) - _helmholtzd(model,T,x)
    f0 = z -> f0β(z[1])
    β =  Optim.optimize(f0,0,βmax).minimizer
    x1 = y
    x2 = (x - β*y)/(1-β)
    X0 = vcat(β,y...,(1-β),(x - β*y)/(1-β)...)
    #println((f0(β),βmax))
        
    
    return (β,x1,x2)
end





#evaluate a function that accepts a vector of length n, N times. this can be ForwardDiffed and accepted by an
#optimization routine
function multi_fval(fn,X,n,N)
res = 0.0
for i = 1:N
    x0 = (i-1)*n+1
    xn = i*n
    res+=fn(X[x0:xn])
end
return res
end

#evaluates helmholtz equations of the form (X,α)
#and adds the term f(X0-X)*(1-sum(α))
#this enforces the restriction of phases and matter conservation
function _multi_fval_red(fn,X,n,N,bufferₓ,X0)
    res = 0.0
    bufferₓ .= 0.0
    α₀=0.0

    for i = 1:N
        x0 = (i-1)*n+1
        xn = i*n
        bufferₓ .= bufferₓ .+ X[x0:xn]
        α₀+=X[end-n+i]
        res+=fn(X[x0:xn])*X[end-n+i]
    end
    return res + fn(X0 .- bufferₓ)*(1-α₀)
    end



