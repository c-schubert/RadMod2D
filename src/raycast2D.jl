#TODO:
# -clean up code (remove unnecessary comments)
# -add docstrings
# -add more tests

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



function get_start_tile(p1::Point2D{T1}, dir::Index2D{T2}, n::T2, dx::T1, 
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



function tilewalk_with_check_for_occ_tiles(fig, ax, p1::Point2D{T1}, p2::Point2D{T1}, 
        dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}
        )::Nothing where {T1<:AbstractFloat, T2<:Integer}

    # do tilewalk between two points and check for occupied tiles
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
    tile = get_start_tile(p1, dir, n, dx, dy)
    println("start tile number: ", tile.x, ", ", tile.y)
    max_steps = get_max_steps(n)
    if !ismissing(t_occ[tile.x, tile.y])
        println("------> occupied ")
        start_walk = false
    else
        # print start tile in color
        poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
        start_walk = true
    end
    if start_walk
        # do walk
        for step = 1:max_steps
            lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
            # check if end reached

            if lx >= (ltot - _TOL) && ly >= (ltot - _TOL)
                break
            end
            tile, lcurrent, pcolor = get_next_tile(tile, dir, lx, ly)
            # println("lcurrent: ",lcurrent)
            # plot next point
            p = p1 + lvec * lcurrent
            scatter!(ax, Point2f(p.x, p.y), color = pcolor, markersize = 15)
            println("next tile number: ", tile.x, ", ", tile.y)
            # check this tile
            if !ismissing(t_occ[tile.x, tile.y])
                println("------> occupied ")
                # plot hitting point
                p = p1 + lvec * lcurrent
                scatter!(ax, Point2f(p.x, p.y), color = :blue, markersize = 15)
                break
            end
            # print tile in color
            poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
        end
    end
end



function tilewalk_with_return!(tile_list::Vector{Index2D{T2}}, p1::Point2D{T1}, 
        p2::Point2D{T1}, dx::T1, dy::T1, n::T2)::T2 where {T1<:AbstractFloat, T2<:Integer}
    # do tilewalk between two points and check for occupied tiles
    max_steps = get_max_steps(n)
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)

    tile = get_start_tile(p1, dir, n, dx, dy)
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