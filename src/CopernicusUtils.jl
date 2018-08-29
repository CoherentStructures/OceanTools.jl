module CopernicusUtils

    using Tensors, StaticArrays, SparseArrays

    import Images
    import NetCDF

    include("interpolation.jl")
    include("dataformat.jl")
    include("misc.jl")


    include("exports.jl")
end


