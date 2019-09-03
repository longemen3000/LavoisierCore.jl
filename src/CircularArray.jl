
struct CircularArray{T,N} <: AbstractArray{T,N}
    x::AbstractArray{T,N}
    function CircularArray(x::AbstractArray{T,N}) where {T,N}
        return new{T,N}(x)
    end
end

Base.size(A::CircularArray) = size(A.x) #size of the vector

Base.length(A::CircularArray)=length(A.x)

function Base.getindex(A::CircularArray, I::Vararg{Int, N}) where N
I2 = size(A)
return Base.getindex(A.x,(mod1.(I,I2))...)
end

Base.getindex(A::CircularArray, I) = [A[i] for i in I]

function Base.setindex!(A::CircularArray,value,I...) 
    I2 = size(A)
    return Base.setindex!(A.x,value,(mod1.(I,I2))...)
end

Base.iterate(S::CircularArray, state=1) = Base.iterate(S.x)
Base.IndexStyle(::Type{CircularArray}) = IndexCartesian()
Base.checkbounds(A::CircularArray, I...) = nothing