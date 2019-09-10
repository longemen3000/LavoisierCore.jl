using LinearAlgebra, JuliaDB, Unitful



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


function molar_to_weight(value,x,MW)
    return value/LinearAlgebra.dot(MW,x)*1000.0
end

function weight_to_molar(value,x,MW)
    return (value*LinearAlgebra.dot(MW,x))/1000.0
end



        
##DIPPR search utils


    
    function property(data::T,name::Symbol) where T <: JuliaDB.Dataset
    return JuliaDB.select(data,name)
    end
    function property(data::T,name::Symbol) where T <: NamedTuple
        return getfield(data,name)
    end
    function property(data::T,name::Symbol) where T <: AbstractDict
        return data[name]
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
 
    """
    loaddippr()
    loads the DIPPR801 database
    # Examples
    ```julia-repl
    julia> molar_to_weight(55.01, [18.01,16.0] [0.5,0.5])
    1
"""
 loaddippr()=JuliaDB.load("DIPPR801.juliadb")



#transformation of arbitrary molar units to the SI units
@inline function _transform_v(v,mw,x)
    return 1.0*v
end


@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐍^-1},mw,x)
    return Unitful.ustrip(u"m^3/mol"(v))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐍^1},mw,x)
    return 1/Unitful.ustrip(u"mol/m^3"(v))
end

@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^-3*Unitful.𝐌^1},mw,x)
    return weight_to_molar(1/Unitful.ustrip(u"kg/m^3"(v)),mw,x)
end


@inline function _transform_v(v::Unitful.Quantity{<:Real, Unitful.𝐋^3*Unitful.𝐌^-1},mw,x)
    return weight_to_molar(Unitful.ustrip(u"m^3/kg"(v)),mw,x)
end

@inline function _transform_T(T)
    return 1.0*T
end

@inline function _transform_T(T::TT) where TT <: Unitful.Temperature
        return Unitful.ustrip(u"K"(T))
end

function _transformVT(V,T,mw,x)
    return (_transform_v(V,mw,x),_transform_T(T))
end

###
### fraction utilities
###
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
function randfrac(N::Integer)
    x = rand(N)
    return x ./sum(x)
end
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
