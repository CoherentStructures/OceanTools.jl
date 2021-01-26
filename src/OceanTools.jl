module OceanTools

    using StaticArrays, SparseArrays, SharedArrays

    import FileIO
    import Tensors
    import Dates
    
    import NCDatasets
    const NCD = NCDatasets

    include("interpolation.jl")
    include("dataformat.jl")
    include("misc.jl")


    include("export.jl")
end
