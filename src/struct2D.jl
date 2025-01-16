################################################################
# Code part of FARADS2D                                        #
# by dominik.bueschgens                                        #
# March 2021                                                   #
################################################################


abstract type Pair2D{T<:Real} <: AbstractVector{T} end

struct Point2D{T<:AbstractFloat} <: Pair2D{T}
    x::T
    y::T
end

struct Vector2D{T<:AbstractFloat} <: Pair2D{T}
    x::T
    y::T
end

struct Index2D{T<:Integer} <: Pair2D{T}
    x::T
    y::T
end


Point2D(x::Real, y::Real) = Point2D(promote(x,y)...)
Point2D(x::Integer, y::Integer) = Point2D(convert(typeof(1.0),x), convert(typeof(1.0),y))

Vector2D(x::Real, y::Real) = Vector2D(promote(x,y)...)
Vector2D(x::Integer, y::Integer) = Vector2D(convert(typeof(1.0),x), convert(typeof(1.0),y))

Index2D(x::Real, y::Real) = Index(convert(typeof(1),x), convert(typeof(1),y))


import Base: convert

convert(::Type{Point2D}, x::Tuple{T,T}) where T <: Number = Point2D(x[1],x[2])
convert(::Type{Vector2D}, x::Tuple{T,T}) where T <: Number = Vector2D(x[1],x[2])
convert(::Type{Index2D}, x::Tuple{T,T}) where T <: Number = Index2D(x[1],x[2])

import Base: size, getindex, unsafe_getindex, setindex!, unsafe_setindex!, +, -, *, /

size(v::Pair2D) = (2, )

# looks like loosing perf if used as vector
getindex(v::Pair2D, i::Int) = i < 3 ? (i == 1 ? v.x : v.y) : error("max index = 2")

unsafe_getindex(v::Pair2D, i::Int) = (i == 1 ? v.x : v.y)

function setindex!(v::Pair2D, a::Real, i::Integer)
    if i < 3
        if i == 1
            v.x = a
        else
            v.y = a
        end
    else
        error("max index = 2")
    end
    return v
end


function unsafe_setindex!(v::Pair2D, a::Real, i::Integer) 
    if i == 1
        v.x = a
    else
        v.y = a
    end

    return v
end

+(a::T, b::T) where T <: Pair2D = T(a.x+b.x,a.y+b.y)
-(a::T, b::T) where T <: Pair2D = T(a.x-b.x,a.y-b.y)
*(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x*b,a.y*b)
/(a::T1, b::T2) where {T1 <: Pair2D, T2 <: Real} = T1(a.x/b,a.y/b)


"""
    norm(v::Point2D)

l-2 norm of vector

#Arguments
- `v::Point2D`: Point2D vector
#Returns
- `l<:AbstractFloat`: length of vector
"""
function norm(v::Point2D{T})::T where T<:AbstractFloat
    # l = LinearAlgebra.norm([v.x, v.y], 2.0)
    l = sqrt(v.x^2 + v.y^2)
    return l
end

function norm(v::Vector2D{T})::T where T<:AbstractFloat
    # l = LinearAlgebra.norm([v.x, v.y], 2.0)
    l = sqrt(v.x^2 + v.y^2)
    return l
end



"""
    normit(v::Point2D)::Point2D

Return normalized vector coordinates as Point2D

#Arguments
- `v::Point2D`: Point2D vector
#Returns
- `vn::Point2D`: With normalized vector coordinates
"""
function normit(v::Point2D{T})::Point2D{T} where T<:AbstractFloat
    vn = v / norm(v)
    return vn
end

function normit(v::Vector2D{T})::Vector2D{T} where T<:AbstractFloat
    vn = v / norm(v)
    return vn
end



"""
    get_com(p1::Point2D{T}, p2::Point2D{T})::Point2D{T} where T<:AbstractFloat

#Arguments
- `p1::Point2D{T}`: first point of line
- `p2::Point2D{T}`: second point of line

#Returns
- `c::Point2D{T}`:  The center of line ((center of mass)) defined by two points p1 and p2 
"""
function get_com(p1::Point2D{T}, p2::Point2D{T})::Point2D{T} where T<:AbstractFloat
    c = p1 + ((p2 - p1) * 0.5)
    return c
end



function get_norm_vec(p1::Point2D{T}, p2::Point2D{T}, dir::Symbol)::Point2D{T} where T<:AbstractFloat
    # calculate normal vector
    n  = Point2D((-1.0) * (p2.y - p1.y), (p2.x - p1.x))
    (dir == :neg) && (n = n * (-1.0))
    n_norm = normit(n)
    return n_norm
end



@inline function get_length(p1::Point2D{T}, p2::Point2D{T})::T where T<:AbstractFloat
    return norm(p1-p2)
end



function get_angle(v1::Point2D{T},v2::Point2D{T})::T where T<:AbstractFloat
    # calculate angle between two vectors
    EPSILON = 1E-8
    cosval::T = dot(v1,v2) / (norm(v1) * norm(v2))
    if cosval < -1 || cosval > 1
        if cosval > (-1 - EPSILON) && cosval < -1
            cosval = -1
        elseif cosval < (1 + EPSILON) && cosval > 1
            cosval = 1
        else
            error("custom error message: cosval in fct: get_angle")
        end
    end
    phir = acos(cosval)
    #phid = phir * 180 / pi
    return phir
end