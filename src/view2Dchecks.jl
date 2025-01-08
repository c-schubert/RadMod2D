"""
    get_number_of_elements_in_tiles(t_occ)

"""
function get_number_of_elements_in_tiles(t_occ)
    nx = size(t_occ,1)
    ny = size(t_occ,2)
    t_num = zeros(Int64, nx, ny)
    for i = 1:nx, j = 1:ny
        if ismissing(t_occ[i,j])
            t_num[i,j] = 0
        else
            t_num[i,j] = size(t_occ[i,j],1)
        end
    end
    return t_num
end



function tile_occ_analysis(t_occ; printit = true)
    # analysis of tile occupation matrix
    t_num = get_number_of_elements_in_tiles(t_occ)
    t_num_v = vec(t_num)
    t_num_v_occ = t_num_v[t_num_v[:] .> 0]
    t_min = minimum(t_num_v_occ)
    t_max = maximum(t_num_v_occ)
    n_all = size(t_num_v,1)
    n_occ = size(t_num_v_occ,1)
    t_mean = sum(t_num_v_occ) / n_all
    t_mean = round(t_mean, digits=1)
    t_mean_occ = sum(t_num_v_occ) / n_occ
    t_mean_occ = round(t_mean_occ, digits=1)
    if printit
        println("Tile occupation with elements:")
        println("    Tiles occupied: ", n_occ, " / ", n_all)
        println("    Max: ", t_max, " // Min: ", t_min)
        println("    Mean: ", t_mean, " // Mean_occ: ", t_mean_occ)
    else
        return t_max, t_min, t_mean_occ, n_occ/n_all
    end 
end