"""
line_segment_intersection(p1::Point2D{T}, p2::Point2D{T}, p3::Point2D{T}, 
        p4::Point2D{T})::Bool where T<:AbstractFloat

Check if two 2D line segments (p1 -> p2) and (p3 -> p4) are intersecting each other.

#Arguments
- `p1::Point2D{T}`: first point of first line segment
- `p2::Point2D{T}`: second point of first line segment
- `p3::Point2D{T}`: first point of second line segment
- `p4::Point2D{T}`: second point of second line segment

#Returns
- `hit::Bool`: true if line segments are intersecting, false otherwise
"""
function line_segment_intersection(p1::Point2D{T}, p2::Point2D{T}, p3::Point2D{T}, 
        p4::Point2D{T})::Bool where T<:AbstractFloat
    # check if two 2D line segments are intersecting
    hit = false
    x1 = p1.x
    x2 = p2.x
    x3 = p3.x
    x4 = p4.x
    y1 = p1.y
    y2 = p2.y
    y3 = p3.y
    y4 = p4.y
    denom = ((y2-y1)*(x3-x4)-(x2-x1)*(y3-y4))
    if denom != 0
        t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/denom
        u = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom
        if (t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0)
            hit = true
        end
    end
    return hit
end


"""
    is_overlapping2D(box1::Vector{T}, box2::Vector{T}) where T <: AbstractFloat

Checks if two 2D bounding boxes are overlapping.

#Arguments
- `box1::Vector{T}`: first bounding box = [xmin, xmax; ymin, ymax]
- `box2::Vector{T}`: second bounding box = [xmin, xmax; ymin, ymax]

"""
# function is_overlapping2D(box1::Vector{T}, box2::Vector{T})::Bool where T <: AbstractFloat
function is_overlapping2D(box1::SVector{4,T}, box2::SVector{4,T})::Bool where T <: AbstractFloat

    # @no_escape begin 
        box1x = @view box1[1:2]
        box1y = @view box1[3:4]
        box2x = @view box2[1:2]
        box2y = @view box2[3:4]

        # box1x = @alloc(T, 2)
        # box1y = @alloc(T, 2)
        # box2x = @alloc(T, 2)
        # box2y = @alloc(T, 2)

        # box1x[1] = box1[1]
        # box1x[2] = box1[2]
        # box1y[1] = box1[3]
        # box1y[2] = box1[4]
        # box2x[1] = box2[1]
        # box2x[2] = box2[2]
        # box2y[1] = box2[3]
        # box2y[2] = box2[4]

        inside = false

        if is_overlapping1D(box1x, box2x) && is_overlapping1D(box1y, box2y)
            inside = true
        end
    # end

    return inside
end



"""
    is_overlapping1D(a1::Vector{T}, a2::Vector{T})::Bool where T <: AbstractFloat

Checks if two 1D bounding "areas" are overlapping. Areas are defined by box1 = [xmin, xmax]
and box2 = [xmin, xmax].

#Arguments
- `box1::Vector{T}`: first bounding line = [xmin, xmax]
- `box2::Vector{T}`: second bounding line = [xmin, xmax]

#Returns
- `inside::Bool`: true if 1D "areas" are overlapping, false otherwise
"""
# function is_overlapping1D(box1::Vector{T},box2::Vector{T})::Bool where T <: AbstractFloat
function is_overlapping1D(box1::SubArray,box2::SubArray)::Bool

    xmin1 = box1[1]
    xmax1 = box1[2]
    xmin2 = box2[1]
    xmax2 = box2[2]

    inside::Bool = false

    if xmax1 >= xmin2 && xmin1 <= xmax1
        inside = true
    end

    return inside
end



"""
    get_area_of_part(m::AbstractModel, idx_part::Integer)

Calculates the area of a part.

#Arguments
- `m::AbstractModel`: model
- `idx_part::Integer`: index of part

#Returns
- `area`: area of part

"""
function get_area_of_part(m::ConstModel, idx_part::Integer)::AbstractFloat
    # returns the area of part p
    start_element_of_part = m.elem2par[idx_part].first
    last_element_of_part = m.elem2par[idx_part].last
    # println("calculate area between elements ", e1, " and ", e2)
    area = 0.0
    for i = start_element_of_part:last_element_of_part
        area += m.elements[i].length
    end

    return area
end



"""
    are_elements_facing(e1::ConstLineElement, e2::ConstLineElement)::Bool

Checks if two elements are facing towards each other. Returns true if the angle between
the normal vector of the elements e1 and e2 and the connector is smaller than 90°.

#Arguments
- `e1::ConstLineElement`: first element
- `e2::ConstLineElement`: second element
"""
function are_elements_facing(e1::ConstLineElement{T1,T2}, e2::ConstLineElement{T1,T2})::Bool where {T1<:Integer, T2<:AbstractFloat} 
    con = e2.com - e1.com
    phir1 = get_angle(e1.norm_vec,con)
    phir2 = get_angle(e2.norm_vec,con*(-1.0))

    if phir1 < (pi * 0.5 - _TOL) && phir2 < (pi * 0.5 - _TOL)
        return true
    end  
    
    return false
end



"""
    existing_vf!(m::ConstModel, vfmat::Matrix{T})::Nothing where T <: AbstractFloat

Checks if element pairs are facing towards each, which means that they theoretically can see
each other if nothing is blocking them and sets the corresponding entry in the passed
blocking matrix to 1 (no blocking), if so. Should be used as initializer for the other
blocking methods.
"""
function existing_vf!(m::ConstModel{T1, T2}, 
            vfmat::Matrix{T2})::Nothing where {T1<:Integer, T2<:AbstractFloat}

    for i1 = 1:m.no_elements, i2 = (i1+1):m.no_elements
        if are_elements_facing(m.elements[i1], m.elements[i2])
            vfmat[i1,i2] = 1.0
            vfmat[i2,i1] = 1.0
        end
    end

    return nothing
end



"""
    blocking_vf_brute_force!(m::ConstModel, vfmat::vfmat{T})::Nothing where T <: AbstractFloat

Checks if element pairs are blocked by others with brute force method. Sets the
corresponding entry in the passed blocking matrix to 0 (blocking), if so. This method is
very slow (O ~ n^3) and should only be used for small models. 
"""
function blocking_vf_brute_force!(m::ConstModel{T1, T2}, 
    vfmat::Matrix{T2})::Nothing where {T1<:Integer, T2<:AbstractFloat}

    @inbounds for i1 = 1:m.no_elements, i2 = (i1+1):m.no_elements
        if isapprox(vfmat[i1,i2], 1.0)
            for ib in 1:m.no_elements
                if ib != i1 && ib != i2
                    p1 = m.elements[i1].com
                    p2 = m.elements[i2].com
                    p3 = m.nodes[m.elements[ib].no_node1]
                    p4 = m.nodes[m.elements[ib].no_node2]
                    blocked = line_segment_intersection(p1, p2, p3, p4)
                    if blocked
                        vfmat[i1,i2] = 0
                        vfmat[i2,i1] = 0
                        break
                    end
                end
            end
        end
    end
end




function blocking_vf_with_tiles_simplified!(m::ConstModel{T2, T1}, vfmat::Matrix{T1}, dx::T1, 
            dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}; elem_in_t = 12::T2, 
            skipped_t::T2 = 2)::Nothing where {T1<:AbstractFloat, T2<:Integer}
    # check if element pairs are blocked by others with tiles
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    count_simp::T2 = 0
    count_all::T2 = 0
    @inbounds for i1 = 1:m.no_elements, i2 = (i1+1):m.no_elements
        # println("--------> checking ",i1," and ",i2)
        p1 = m.elements[i1].com
        p2 = m.elements[i2].com
        if isapprox(vfmat[i1,i2],1.0)
            # first tilewalk and get all tiles between i1 and i2
            ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
            hitten = false
            @inbounds for i = 1:ntiles
                t = tiles[i]
                # println("next tile number: ", t.x, ", ", t.y)
                if !ismissing(t_occ[t.x,t.y])
                    # println("------> occupied ")
                    count_all += 1
                    if (size(t_occ[t.x,t.y],1) <= elem_in_t || i <= 1+skipped_t || 
                                i >= ntiles-skipped_t)

                        for it in t_occ[t.x,t.y]
                            if it != i1 && it != i2
                                # println("------> ",i1," and ",i2, " checking with ",it)
                                p3 = m.nodes[m.elements[it].no_node1]
                                p4 = m.nodes[m.elements[it].no_node2]
                                blocked = line_segment_intersection(p1, p2, p3, p4)
                                if blocked
                                    # println("----------------> ",i1," and ",i2," intersected by ",it)
                                    vfmat[i1,i2] = 0.0
                                    vfmat[i2,i1] = 0.0
                                    hitten = true
                                    break
                                end
                            end
                        end
                    else
                        # println("----------------> SIMPLIFIED -> ",i1," and ",i2," intersected by ",it)
                        vfmat[i1,i2] = 0.0
                        vfmat[i2,i1] = 0.0
                        count_simp += 1
                    end
                end
                if hitten
                    break
                end
            end
        end
    end
    # println("counter simp/all: ", count_simp, " / ", count_all)
    println("simplification rate: ", round(count_simp/count_all*100, digits = 2), " %")

    return nothing
end



"""
    calculating_vf!(m::ConstModel, mat::Matrix{Int64}; normit = false)

Calculates the viewfactors between all element combunatations i-j where mat[j][j] == 1
in the model `m` and stores them in the matrix `mat`.

# Arguments
- `m::ConstModel`: model
- `mat::Matrix{Int64}`: matrix with viewfactors

# Keyword Arguments
- `normit::Bool = false`: if true, the viewfactors are normalized to 1

# Returns
- `nothing`: the viewfactors are updated within the matrix `mat`
"""
function calculating_vf!(m::ConstModel{T2, T1}, mat::Matrix{T1}; normit::Bool = false
            )::Nothing  where {T1<:AbstractFloat, T2<:Integer}
    # calculate viewfactors if vf is existing
    for i1 = 1:m.no_elements, i2 = (i1+1):m.no_elements
        if isapprox(mat[i1,i2],1.0)
            # viewfactor existing
            a = m.nodes[m.elements[i1].no_node1]
            b = m.nodes[m.elements[i1].no_node2]
            c = m.nodes[m.elements[i2].no_node1]
            d = m.nodes[m.elements[i2].no_node2]
            # punkte müssen richtig sortiert sein
            # nur aus python rauskopiert - nicht mehr kontrolliert
            # hier umgangen mit betrag
            # 4 Strecken berechnen -> Abstand zwischen 2 Punkten
            ac = norm(c - a)
            bd = norm(d - b)
            bc = norm(c - b)
            ad = norm(d - a)
            # println("AC: ", ac," BD: ", bd," BC: ", bc," AD: ", ad)
            # cross strings formula - Fläche berechnen und vf bestimmen
            A1 = norm(b - a)
            A2 = norm(d - c)
            mat[i1,i2] = abs((0.5 * ((ac+bd)-(bc+ad))) / A1)
            mat[i2,i1] = abs((0.5 * ((ac+bd)-(bc+ad))) / A2)
            #print(i1," -> ",i2,":",vf12,vf21)
        end
    end
    if normit
        control = sum(mat, dims=2)
        mat .= mat ./ control
    end

    return nothing
end


