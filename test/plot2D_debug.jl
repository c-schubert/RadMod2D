
function blocking_vf_with_tiles_2elem(fig, ax, i1, i2, m::ConstModel, dx, dy, n::T2, t_occ) where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elements[i1].com
    p2 = m.elements[i2].com
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    # first tilewalk and get all tiles between i1 and i2
    ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
    println(tiles[1:ntiles])
    hitten = false
    for i = 1:ntiles
        t = tiles[i]
        println("next tile number: ", t.x, ", ", t.y)
        # print tile in color
        poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:grey, 0.6))
        if !ismissing(t_occ[t.x,t.y])
            println("------> occupied ")
            for it in t_occ[t.x,t.y]
                if it != i1 && it != i2
                    println("------> ",i1," and ",i2, " checking with ",it)
                    p3 = m.nodes[m.elements[it].no_node1]
                    p4 = m.nodes[m.elements[it].no_node2]
                    blocked = line_segment_intersection(p1, p2, p3, p4)
                    if blocked
                        println("----------------> ",i1," and ",i2," intersected by ",it)
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

function blocking_vf_with_tiles_2elem_tiles(fig, ax, i1, i2, m::ConstModel, dx, dy, n::T2, t_occ)::Nothing where T2<:Integer
    # check if element pair is blocked by others with tiles
    p1 = m.elements[i1].com
    p2 = m.elements[i2].com
    max_steps = get_max_steps(n)
    tiles = Vector{Index2D{T2}}(undef,max_steps)
    # first tilewalk and get all tiles between i1 and i2
    ntiles = tilewalk_with_return!(tiles, p1, p2, dx, dy, n)
    # println(tiles[1:ntiles])
    for i = 1:ntiles
        t = tiles[i]
        println("next tile number: ", t.x, ", ", t.y)
        # print tile in color
        poly!(ax, Rect((t.x-1) * dx, (t.y-1) * dy, dx, dy), color = (:grey, 0.6))
    end

    return nothing
end


function tilewalk_with_return(fig, ax, p1::Point2D{T1}, p2::Point2D{T1}, dx::T1, dy::T1, 
    n::T2)::Vector{Index2D{T2}} where {T1<:AbstractFloat, T2<:Integer}
# do tilewalk between two points and check for occupied tiles
    max_steps = get_max_steps(n)
    tile_list = Vector{Index2D{T2}}(undef,max_steps)
    dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
    tile = find_grid_point(p1, dir, n, dx, dy)
    println("start tile number: ", tile.x, ", ", tile.y)
    hit = 1
    tile_list[hit] = Index2D(tile.x, tile.y)
    # print start tile in color
    poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
    # do walk
    for step = 1:max_steps
        lx, ly = get_lengths_to_next_tiles(p1, p2, dir, tile, dx, dy)
        # check if end reached

        if lx >= (ltot - _TOL) && ly >= (ltot - _TOL)
            break
        end
        tile, lcurrent, pcolor = get_next_tile(tile, dir, lx, ly)
        println("next tile number: ", tile.x, ", ", tile.y)
        hit += 1
        tile_list[hit] = Index2D(tile.x, tile.y)
        # print tile in color
        # just for graphical debuging moved out of package
        poly!(ax, Rect((tile.x-1) * dx, (tile.y-1) * dy, dx, dy), color = (:blue, 0.3))
    end
    return tile_list[1:hit]
end



function tilewalk_with_check_for_occ_tiles(fig, ax, p1::Point2D{T1}, p2::Point2D{T1}, 
    dx::T1, dy::T1, n::T2, t_occ::Matrix{Union{Vector{T2},Missing}}
    )::Nothing where {T1<:AbstractFloat, T2<:Integer}

# do tilewalk between two points and check for occupied tiles
dir, ltot, lvec = get_vectors_for_tilewalk(p1, p2)
tile = find_grid_point(p1, dir, n, dx, dy)
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

