module CopernicusUtils

    using Tensors, StaticArrays, SparseArrays, SharedArrays

    import Images
    import NCDatasets
    import NetCDF
    const NCD = NCDatasets

    include("interpolation.jl")
    include("dataformat.jl")
    include("misc.jl")


    include("export.jl")
end
