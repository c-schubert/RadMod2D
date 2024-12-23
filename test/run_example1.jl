using RadMod2D

# model
include("./models2D.jl")
r1x = 0.8
r1y = 0.8
r2x = 1.2
r2y = 1.2
elemsize = 0.05
m = model_rect_in_rect(r1x, r1y, r2x, r2y, elemsize)

# view factors
vfmat = zeros(Float64, m.no_elements, m.no_elements)
@time existing_vf!(m, vfmat)
n = 15
dx, dy = get_tile_dimensions(m, n)
@time t_occ = check_tile_occupation(m, dx, dy, n)
@time blocking_vf_with_tiles!(m, vfmat, dx, dy, n, t_occ)
@time calculating_vf!(m, vfmat, normit = true)
vfmatp = compact_vfmat_to_parts(m, vfmat, normit = true)
println("vf mat for parts:")
for i = 1:size(vfmatp,1)
    println(vfmatp[i,:])
end

# net radiation method
epsilon = zeros(m.no_elements,1)
set_bc_part!(m, epsilon, 1, 0.3)
set_bc_part!(m, epsilon, 2, 0.6)
temp = zeros(m.no_elements,1)
set_bc_part!(m, temp, 1:2, 300)
temp[79:106,1] .= 600 # bc for specific elements
@time Q, G = tempsolver(m, vfmat, temp, epsilon)
length = [m.elements[i].length for i = m.elem2par[1].first:m.elem2par[end].last]
q = Q[:] ./ length[:]

println("calculation finished")

# include("./test/run_example1.jl")