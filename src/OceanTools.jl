module OceanTools

    using Tensors, StaticArrays, SparseArrays, SharedArrays

    import FileIO
    import NCDatasets
    import Dates
    const NCD = NCDatasets

    include("interpolation.jl")
    include("dataformat.jl")
    include("misc.jl")


    include("export.jl")
end
