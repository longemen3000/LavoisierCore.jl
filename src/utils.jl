using LinearAlgebra, JuliaDB



"""
    disallowmissing(x)
    transforms an array of type Array{Union{Missing,T},N}, to Array{T,N},
    throws an error otherwise
    # Examples
    ```julia-repl

""" 
function disallowmissing(x::Array{Union{Missing,T},N}) where T where N
    if any(ismissing,x)
        throw(DomainError(x, """
        
        a missing was detected. this model requires 
        that all the values of the input properties should exist.

        """))
    else
        return convert(Array{T,N},x)
    end
end

"""
    molar_to_weight(value,x,MW)
    transforms from 'x/mol' to 'x/kg', using the molar fraction 'x' 
    and the molecular weight 'MW' 
    # Examples
    ```julia-repl
    julia> molar_to_weight(55.01, [18.01,16.0] [0.5,0.5])
    1
""" 
function molar_to_weight(value,x::Array{T,1},MW::Array{T,1}) where T <: Real
    return value/LinearAlgebra.dot(MW,x)/1000
end
"""
    weight_to_molar(value,x,MW)
    transforms from 'x/kg' to 'x/mol', using the molar fraction 'x' 
    and the molecular weight 'MW' 
    # Examples
    ```julia-repl
    julia> molar_to_weight(55.01, [18.01,16.0] [0.5,0.5])
    1
"""
function weight_to_molar(value::T,x::Array{T,1},MW::Array{T,1}) where T <: Real
    return value*LinearAlgebra.dot(MW,x)*1000
end

@inline function reducedvolume(x,rhoc,betav,gammav)
N = length(x)
x2 = copy(x)

res1 = x[N]^2/rhoc[N]
for i = 1:N-1
    
    res1 += x[i]^2/rhoc[i]
    for j=(i+1):N
        res1 += 0.25*x[i]*x[j]*
        betav[i,j]*
        gammav[i,j]*
        ((x[i]+x[j])/(x[i]*betav[i,j]^2 + x[j])) *
        (rhoc[i]^(-1/3)+rhoc[j]^(-1/3))
    end
end

return 1.0/res1
end

@inline function reducedvolume(x,rhoc)
    N = length(x)
    res1 = x[N]^2/rhoc[N]
    for i = 1:N-1
        res1 += x[i]^2/rhoc[i]
         for j=(i+1):N
            res1 += 0.25*x[i]*x[j]*
            (rhoc[i]^(-1/3)+rhoc[j]^(-1/3))
        end
    end
    
    return 1.0/res1
end

@inline function reducedtemperature(x,Tc,betat,gammat)
N = length(x)
res1 = x[N]^2*Tc[N]
for i = 1:N-1
    res1 += x[i]^2*Tc[i]
    for j=(i+1):N
        res1 += 0.25*x[i]*x[j]*
        betat[i,j]*
        gammat[i,j]*
        ((x[i]+x[j])/(x[i]*betat[i,j]^2 + x[j])) *
        sqrt(Tc[i]*Tc[j])
    end
end
return res1
end

@inline function reducedtemperature(x,Tc)
    N = length(x)
    res1 = x[N]^2*Tc[N]
    for i = 1:N-1
        res1 += x[i]^2*Tc[i]
        for j=(i+1):N
            res1 += 0.25*x[i]*x[j]*sqrt(Tc[i]*Tc[j])
        end
    end
    return res1
end

function normalizefrac!(x)
    summ = sum(x)
    for i = 1:length(x)
        x[i] /=summ
    end
end

function normalizefrac(x)
    summ = sum(x)
    y = copy(x)
    for i = 1:length(x)
        y[i] /=summ
    end
    return y
end
        
##DIPPR search utils

function ischemicaldata(data::T) where T <: JuliaDB.Dataset
    return true
end
    
    function property(data::T,name::Symbol) where T <: JuliaDB.Dataset
    return JuliaDB.select(data,name)
    end
    
    function compound(db::T,index::Int) where T <: JuliaDB.Dataset
        return db[index]
    end
    

    function searchdata(db::T,value::String) where T <: JuliaDB.Dataset
    names = (:name,:synonyms,:iupacName,:casName,:synonyms)
    indxs = Array{Int64,1}()
        for i = 1:5 
        x = JuliaDB.rows(db, names[i])
        indxs = vcat(findall(z->occursin(value,z),x),indxs)
        end
        indxs=unique!(indxs)
    return (JuliaDB.subtable(db, indxs),indxs)
    end
    
 