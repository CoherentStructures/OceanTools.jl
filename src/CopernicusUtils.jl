module CopernicusUtils

    using Tensors, StaticArrays, SparseArrays, SharedArrays

    import Images
    import NetCDF

    include("interpolation.jl")
    include("dataformat.jl")
    include("misc.jl")


    include("export.jl")
end
