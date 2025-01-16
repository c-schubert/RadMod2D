using RadMod2D
using Statistics

include("./plot2D.jl")
include("./plot2D_debug.jl")

function punkt_auf_kreis(r::Float64, versatz_grad::Float64)
    # Umwandlung von Grad in Bogenmaß
    versatz_rad = deg2rad(versatz_grad)
    
    # Berechnung der x- und y-Koordinaten
    x = r * sin(versatz_rad)
    y = r * cos(versatz_rad)
    
    return (x, y)
end

tube_diameter = 0.168
heater_diameter = 0.007
inner_rod_diameter = 0.03

inner_rad = 0.039
outer_rad = 0.05975

no_outer = 18
no_inner = 12

deg_start_center_outer = 0
deg_start_center_inner = 5

inner_center = Point2D[]
outer_center = Point2D[]

println(deg_start_center_inner)

m = create_empty_model()
center = Point2D(0.0,0.0)

deg_offset_outer = 360/(no_outer)
deg_offset_inner = 360/(no_inner)

for i = 1:no_outer
    push!(outer_center, punkt_auf_kreis(outer_rad, deg_offset_outer * (i-1) + deg_start_center_outer))
end

for i = 1:no_inner
    push!(inner_center, punkt_auf_kreis(inner_rad, deg_offset_inner * (i-1) + deg_start_center_inner))
end

n_lines_heater =round(Int64,70/2)
n_lines_tube = round(Int64,1680/4)
n_lines_rod = round(Int64,300/3)

println("E Size heater $(pi*heater_diameter/n_lines_heater)")
println("E Size tube $(pi*tube_diameter/n_lines_tube)")
println("E Size rod $(pi*inner_rod_diameter/n_lines_rod)")

add!(m, circle(tube_diameter, center, seed = n_lines_tube, 
        dir = :pos, name = "Tube"))

add!(m, circle(inner_rod_diameter, center, seed = n_lines_rod, 
        dir = :neg, name = "Rod"))

for i = 1:no_outer
    add!(m, circle(heater_diameter, outer_center[i], seed = n_lines_heater, 
            dir = :neg, name = "Outer Heater $i"))
end
    
for i = 1:no_inner
    add!(m, circle(heater_diameter, inner_center[i], seed = n_lines_heater, 
            dir = :neg, name = "Inner Heater $i"))
end

m_const = ConstModel(m)

fig, ax = new_figure()
plot_model(fig, ax, m_const, show_norm_vec = true, show_nodes = false, show_com = false, norm_vec_scale = 0.0075)
# save("figure_tubothal.png", fig, pdf_version="1.4")


n = 20
dx, dy = get_tile_deltas(m_const, n)
vfmat = zeros(Float64, m_const.no_elements, m_const.no_elements)
existing_vf!(m_const, vfmat)
tile_orgin = get_tile_grid_origin(m_const)

tiles = TileGrid(tile_orgin, Vector2D(dx,dy), Index2D(n,n))

t_occ = check_tile_occupation(m_const, dx, dy, n)


# Plots.spy(ismissing.(t_occ))

# this is working fine
# blocking_vf_brute_force!(m_const, vfmat)

# this is not -> seems to be that tiles does not support negativ model coordinates -> modify tile model!!
blocking_vf_with_tiles!(m_const, vfmat, dx, dy, n, t_occ) # here is a bug somewhere / even if RadMod2D tests are running fine ...

calculating_vf!(m_const, vfmat, normit = false)
vfmatp = compact_vfmat_to_parts(m_const, vfmat, normit = true)

vf_tube_to_tube = vfmatp[1,1]
vf_tube_to_rod = vfmatp[1,2]

vf_tube_to_outer_heaters = sum(vfmatp[1,3:2+no_outer])
vf_tube_to_inner_heaters = sum(vfmatp[1,end-(no_inner-1):end])

avg_vf_outer_heater_to_tube = mean(vfmatp[3:2+no_outer,1])
avg_vf_inner_heaters_to_tube = mean(vfmatp[end-(no_inner-1):end,1])


println("vf_tube_to_tube $vf_tube_to_tube")
println("vf_tube_to_rod $vf_tube_to_rod")
println("vf_tube_to_outer_heaters $vf_tube_to_outer_heaters")
println("vf_tube_to_inner_heaters $vf_tube_to_inner_heaters")
println("avg_vf_outer_heater_to_tube $avg_vf_outer_heater_to_tube")
println("avg_vf_inner_heaters_to_tube $avg_vf_inner_heaters_to_tube")

# save_vfmat_part_to_german_csv(vfmatp, filename="vfmatp_tubo.csv", partmatrix=true, model=m_const)

