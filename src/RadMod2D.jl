module RadMod2D

    # julia packages
    using DelimitedFiles
    using StaticArrays
    using LinearAlgebra
    using FastPow
    using Printf

    
    # TODO: image export functionality within plotting functions
    # using CairoMakie # for image export
    # CairoMakie.activate!(type = "svg")


    const _TOL = 1E-12
    
    include("./struct2D.jl")
    export Point2D, Index2D
    export norm, get_length

    include("./mesh2D.jl")
    export AbstractModel, MutableModel, ConstModel
    export create_empty_model
    export add!, offset_model!
    export ConstModel
    export element_analysis

    include("./geom2D.jl")
    export edge, rectangle, circle, circle_open, cosinus

    include("./raycast2D.jl")
    export create_randomly_occupied_tiles
    export get_max_steps
    export tilewalk_with_check_for_occ_tiles
    export tilewalk_with_return, tilewalk_with_return!

    include("./view2D.jl")
    export existing_vf!
    export blocking_vf_brute_force!
    export blocking_vf_with_tiles!
    export blocking_vf_with_tiles_simplified!
    export calculating_vf!
    export get_area_of_part
    export get_tile_dimensions
    export are_elements_facing
    export check_tile_occupation
    export line_segment_intersection

    include("./view2Dchecks.jl")
    # comapre elements checks
    export tile_occ_analysis
    export get_number_of_elements_in_tiles

    include("./view2Dpostproc.jl")
    export compact_vfmat_to_parts
    export save_vfmat_part_to_csv
    export save_vfmat_part_to_german_csv

    include("./therm2D.jl")
    export sigma
    export set_bc_part!
    export tempsolver

end