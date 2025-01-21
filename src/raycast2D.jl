#TODO:
# -clean up code (remove unnecessary comments)
# -add docstrings
# -add more tests


struct TileGrid{T1<:AbstractFloat, T2<:Integer}
    origin::Point2D{T1}
    del::Vector2D{T1}
    n::Index2D{T2}
end

struct OccupiedTileGrid{T1<:AbstractFloat, T2<:Integer}
    tg::TileGrid{T1, T2}
    occ::Matrix{Union{Vector{T2},Missing}}
end


"""
    get_tile_deltas(m::ConstModel, nx::Integer, ny::Integer)::Tuple

Gets a tiling of the model with nx x ny tiles. Returns the dimensions of a single
tile delta_x and delta_y.
"""
function get_tile_deltas(m::ConstModel{T1, T2}, nx::Integer, ny::Integer)::Tuple{T2,T2} where {T1<:Integer, T2<:AbstractFloat}
    # create tiles based on model nodes
    # get min and max x and y of nodes
    nmin, nmax = get_min_max_coordinates(m)
    delta_x::T2 = (nmax.x - nmin.x) / nx
    delta_y::T2 = (nmax.y - nmin.y) / ny
    return delta_x, delta_y
end

function get_tile_deltas(m::ConstModel{T1, T2}, n::Integer)::Tuple{T2,T2} where {T1<:Integer, T2<:AbstractFloat}
    return get_tile_deltas(m, n, n)
end


"""
    get_tile_grid_origin(m::ConstModel{T1, T2})::Point2D{T2} where {T1<:Integer, T2<:AbstractFloat}

Sets min.x min.y coordinates of the model as tile origin (cartesian coords)
"""
function get_tile_grid_origin(m::ConstModel{T1, T2})::Point2D{T2} where {T1<:Integer, T2<:AbstractFloat}
    nmin, nmax = get_min_max_coordinates(m)

    return Point2D(nmin.x, nmin.y)
end


function get_tile_grid(m::ConstModel{T1, T2}, nx::T1, ny::T1)::TileGrid{T2,T1} where {T1<:Integer, T2<:AbstractFloat}
    origin = get_tile_grid_origin(m)
    delx, dely = get_tile_deltas(m, nx, ny)

    return TileGrid(origin, Vector2D(delx, dely), Index2D(nx, ny))
end

function get_tile_grid(m::ConstModel{T1, T2}, n::T1)::TileGrid{T2,T1} where {T1<:Integer, T2<:AbstractFloat}
    return get_tile_grid(m,n,n)
end


"""
    check_tile_occupation(m::ConstModel{T1, T2}, dx::T2, dy::T2, n::T1; blockparts::Vector{T1} = 1:m.no_parts)::Matrix{Union{Vector{T1},Missing}} where {T1<:Integer, T2<:AbstractFloat}

Returns a matrix of size n x n type Union{Vector, Missing}. If not missing the vector
contains the indices of the elements that are within the tile. 
"""
function check_tile_occupation(m::ConstModel{T1, T2}, tg::TileGrid{T2, T1}; 
            blockparts=1:m.no_parts)::Matrix{Union{Vector{T1},Missing}} where {T1<:Integer, T2<:AbstractFloat}
    
    dx = tg.del.x
    dy = tg.del.y

    nx = tg.n.x
    ny = tg.n.y

    # contains the elements that are in the tile for each tile in n x n tiles
    t_occ = Matrix{Union{Vector{T1}, Missing}}(missing, nx, ny)

    # contains the tile that has the element
    element_inside_tile_no = Vector{T1}(undef, m.no_elements)

    for t1 = 1:nx
        for t2 = 1:ny
            tmin_x = dx*(t1-1) + tg.origin.x
            tmax_x = dx*t1  + tg.origin.x
            tmin_y = dy*(t2-1)  + tg.origin.y
            tmax_y = dy*t2  + tg.origin.y
            hit = 0
            for p in blockparts
                e1 = m.elem2par[p].first
                e2 = m.elem2par[p].last
                for i = e1:e2
                    no_node1 = m.nodes[m.elements[i].no_node1]
                    no_node2 = m.nodes[m.elements[i].no_node2]
                    if is_line_segment_inside_tile_BB(no_node1, no_node2, tmin_x, 
                                tmax_x, tmin_y, tmax_y)
                        hit += 1
                        element_inside_tile_no[hit] = i
                    end
                end
            end
            if hit > 0
                t_occ[t1,t2] = element_inside_tile_no[1:hit]
            end
        end
    end
    return t_occ
end

"""

legacy workaround before TileGrid struct
"""
function check_tile_occupation(m::ConstModel{T1, T2}, dx::T2, dy::T2, n::T1; 
    blockparts=1:m.no_parts)::Matrix{Union{Missing, Vector{T1}}} where {T1<:Integer, T2<:AbstractFloat}

    tile_orgin = get_tile_grid_origin(m)
    tiles = TileGrid(tile_orgin, Vector2D(dx,dy), Index2D(n,n))

    return check_tile_occupation(m, tiles, blockparts=blockparts)
end

import Base: ismissing

function ismissing(tile_occ_matrix::Matrix{Union{Missing, Vector{Integer}}})::Matrix{Bool, Bool}

    missing_mat = falses(size(tile_occ_matrix))

    for i in eachindex(tile_occ_matrix,1)
        for j in eachindex(tile_occ_matrix,2)
            if ismissing(tile_occ_matrix[i,j])
                missing_mat[i,j] = true
            end
        end
    end

    return missing_mat
end


function blocking_vf_with_tiles!(m::ConstModel{T2, T1}, mat::Matrix{T1}, 
        tg::OccupiedTileGrid{T2, T1})::Nothing where {T1<:AbstractFloat, T2<:Integer}

    # check if element pairs are blocked by others with tiles
    max_steps = get_max_steps(tg.n.x, tg.n.y)
    tiles = Vector{Index2D{T2}}(undef,max_steps)

    @inbounds for i1 = 1:m.no_elements, i2 = (i1+1):m.no_elements
        # println("--------> checking ",i1," and ",i2) 
        p1 = m.elements[i1].com
        p2 = m.elements[i2].com
        if isapprox(mat[i1,i2], 1.0)
            # first tilewalk and get all tiles between i1 and i2
            ntiles = tilewalk_with_return!(tiles, p1, p2, tg.del.x, tg.del.y, tg.nx, tg.ny)
            hitten = false
            for i = 1:ntiles
                t = tiles[i]
                println("next tile number: ", t.x, ", ", t.y)
                if !ismissing(t_occ[t.x,t.y])
                    # println("------> occupied ")
                    @inbounds for it in t_occ[t.x,t.y]
                        if it != i1 && it != i2
                            # println("------> ",i1," and ",i2, " checking with ",it)
                            p3 = m.nodes[m.elements[it].no_node1]
                            p4 = m.nodes[m.elements[it].no_node2]
                            blocked = line_segment_intersection(p1, p2, p3, p4)
                            if blocked
                                # println("----------------> ",i1," and ",i2," intersected by ",it)
                                mat[i1,i2] = 0
                                mat[i2,i1] = 0
                                hitten = true
                                break
                            end
                        end
                    end
                end
                if hitten
                    break
                end
            end
        end
    end

    return nothing
end


"""
    get_tilewalk(p1::Point2D{T1}, p2::Point2D{T1}, n::T2, dx::T1, dy::T1)::Vector{Index2D{T2}}
    where {T1<:AbstractFloat, T2<:Integer}

Calculates the tile walk from coordinate points p1 to p2. The tile walk is a vector of
tile indices. The tile walk is calculated by the Bresenham algorithm. Returns the index
for the next tilewalk step. As (x_ind, y_ind) where x_ind and y_ind are the indices for
the next tile in x and y direction. 
"""
function get_vectors_for_tilewalk(p1::Point2D{T1}, 
            p2::Point2D{T1})::Tuple{Index2D{Int64}, T1, Point2D{T1}} where T1<:AbstractFloat
    # calculate necessary data for tilewalk
    vec = p2 - p1
    vec_l = norm(vec)
    vec_norm = normit(vec)

    dirx::typeof(1) = 0
    if vec_norm.x > _TOL
        dirx = 1
    elseif vec_norm.x < ((-1) * _TOL)
        dirx = -1
    end

    diry::typeof(1) = 0
    if vec_norm.y > _TOL
        diry = 1
    elseif vec_norm.y < ((-1) * _TOL)
        diry = -1
    end

    return Index2D(dirx, diry), vec_l, vec_norm
end

function find_grid_point(p1::Point2D{T1}, tg::TileGrid{T2, T1};epsilon = _TOL) where {T1<:AbstractFloat, T2 <: Integer}
    px = p1.x
    py = p2.x

    x0 = tg.origin.x
    y0 = tg.origin.y

    dx = tg.del.x
    dy = tg.del.y

    # Berechnung der Zellenindizes (i, j)
    i = (px - x0) / dx
    j = (py - y0) / dy
    
    # Toleranzprüfung und Runden der Indizes, wenn der Punkt nahe der Gitterposition ist
    if abs(px - (x0 + i * dx)) < epsilon
        i = round(i)  # Runden auf den nächsten Zellindex
    end
    
    if abs(py - (y0 + j * dy)) < epsilon
        j = round(j)  # Runden auf den nächsten Zellindex
    end
    
    # Berechnung der tatsächlichen Position im Gitter
    x = x0 + i * dx
    y = y0 + j * dy
    
    return (x, y), (i, j)
end


function get_start_tile(p1::Point2D{T1}, dir::Index2D{T2}, nx::T2, ny::T2, dx::T1, 
        dy::T1)::Index2D{T2} where {T1<:AbstractFloat, T2 <: Integer}
    # get tile of starting point
    tx_step = 1
    ty_step = 1

    if dir.x > 0
        tx_step += floor(T2, (p1.x + dir.x * _TOL) / dx)
    elseif dir.x < 0
        tx_step += floor(T2, (p1.x + dir.x * _TOL) / dx)
    else
        # wenn exakt auf einer tile Linie -> 1D walk
        # entscheidet Vorzeichen vor _TOL wo der Bucket liegt
        # wichtig wenn linie am Rand
        if p1.x > (dx * n - _TOL)
            tx_step += floor(T2, (p1.x - _TOL) / dx)
        else
            tx_step += floor(T2, (p1.x + _TOL) / dx)
        end
    end
    if dir.y > 0
        ty_step += floor(T2, (p1.y + dir.y * _TOL) / dy)
    elseif dir.y < 0
        ty_step += floor(T2, (p1.y + dir.y * _TOL) / dy)
    else
        # wenn exakt auf einer tile Linie -> 1D walk
        # entscheidet Vorzeichen vor _TOL wo der Bucket liegt
        # wichtig wenn linie am Rand
        if p1.y > (dy * n - _TOL)
            ty_step += floor(T2, (p1.y - _TOL) / dy)
        else
            ty_step += floor(T2, (p1.y + _TOL) / dy)
        end
    end

    # @assert (tx_step > 0 && tx_step <= n)
    # @assert (ty_step > 0 && ty_step <= n)

    return Index2D(tx_step, ty_step)
end



@inline function Sx(m::T)::T where T<:AbstractFloat
    @fastpow val = sqrt(1 + m^2)
    return val
end



@inline function Sy(m::T)::T where T<:AbstractFloat
    @fastpow val = sqrt(1 + (m^(-1))^2)
    return val
end



function get_lengths_to_next_tiles(p1::Point2D{T1}, p2::Point2D{T1}, dir::Index2D{T2}, 
        tile::Index2D{T2}, dx::T1, dy::T1)::Tuple{T1,T1} where {T1<:AbstractFloat, T2<:Integer}
    # get lengths in x and y direction to next tiles
    
    lx::T1 = Inf
    ly::T1 = Inf

    dx_step::T1 = NaN
    dy_step::T1 = NaN

    if dir.x > 0
        dx_step = dx * (tile.x)
    elseif dir.x < 0
        dx_step = dx * (tile.x - 1)
    end

    if dir.y > 0
        dy_step = dy * (tile.y)
    elseif dir.y < 0
        dy_step = dy * (tile.y - 1)
    end

    # correction for 1D case and changing NaN into Inf
    if !isnan(dx_step)
        lx = dir.x * (dx_step - p1.x) * Sx((p2.y-p1.y) / (p2.x-p1.x))
    end

    if !isnan(dy_step)
        ly = dir.y * (dy_step - p1.y) * Sy((p2.y-p1.y) / (p2.x-p1.x))
    end

    # println("lx: ",lx,"  ly: ",ly)
    return lx, ly
end



function get_next_tile(tile::Index2D{T2}, dir::Index2D{T2}, lx::T1, 
        ly::T1)::Tuple{Index2D{T2},T1,Symbol} where {T1<:AbstractFloat, T2<:Integer}
    # get next tile
    # bei einem 2D Fall exakt durch die Knoten entscheidet 
    # hier < und <= welche Seite das next tile liegt
    # default von mir: < geht immer in nächste y richtung

    # TODO Changes for non quad grid?!!!


    if abs(lx) < abs(ly)
        # println("moving in x direction")
        tile_new = tile + Index2D(dir.x, 0)
        lcurrent = lx
        color = :cyan
    else
        # println("moving in y direction")
        tile_new = tile + Index2D(0, dir.y)
        lcurrent = ly
        color = :coral
    end
    return tile_new, lcurrent, color
end


"""
    create_randomly_occupied_tiles(n::T; density = 0.1)::Matrix{Union{Vector{T},Missing}} where T<:Integer

Debug/Test funciton which randomly creates occupied tiles.
"""
function create_randomly_occupied_tiles(n::T; density = 0.1)::Matrix{Union{Vector{T},Missing}} where T<:Integer
    # create random occ tiles
    # only for testing
    t_occ = rand(n,n)
    t_occ[t_occ .> (1-density)] .= 1
    t_occ[t_occ .<= (1-density)] .= 0
    t_occ = convert.(Integer,t_occ)
    t_occ_n = Matrix{Union{Vector{T},Missing}}(missing,n,n)
    # t_occ_n[t_occ .== 1] .= Vector{T}(undef,2)
    # t_occ = rand([0, 1], n, n)
    # t_occ = zeros(Integer, n, n)
    for i = 1:n, j = 1:n
        if isapprox(t_occ[i,j], 1.0)
            t_occ_n[i,j] = [1,1]
        end
    end
    return t_occ_n
end


@inline function get_max_steps(n::T)::T where T<:Integer
    @fastpow m1 = round(T,sqrt(n^2+n^2)*2)
    # m2 = convert(typeof(n),m1)
    return m1
end

@inline function get_max_steps(nx::T,ny::T)::T where T<:Integer
    @fastpow m1 = round(T,sqrt(nx^2+ny^2)*2)
    # m2 = convert(typeof(n),m1)
    return m1
end


function tilewalk_with_return!(tile_list::Vector{Index2D{T2}}, p1::Point2D{T1}, 
        p2::Point2D{T1}, dx::T1, dy::T1, nx::T2, ny::T2)::T2 where {T1<:AbstractFloat, T2<:Integer}
    # do tilewalk between two points and check for occupied tiles
    max_steps = get_max_steps(n)
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)

    tile = get_start_tile(p1, dir, nx, ny, dx, dy)
    hit = 1
    tile_list[hit] = tile
    for step = 1:max_steps # do walk
        if step == max_steps
            @warn "reached max allowed tiles in tilewalk"
        end
        lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
        # check if end reached

        if lx >= (ltot - _TOL) && ly >= (ltot - _TOL)
            break
        end

        tile, nothing, nothing = get_next_tile(tile, dir, lx, ly)
        hit += 1
        
        @assert (tile.x > 0 && tile.x <= n)
        @assert (tile.y > 0 && tile.y <= n)

        tile_list[hit] = tile
    end
    return hit
end