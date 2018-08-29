function earthDiffTensor(x::Vec{2})
    #const R = 6371e3 # Radius of the Earth in metres
    return SymmetricTensor{2,2,Float64,3}((1. /cos(deg2rad(x[2]))^2,0.0,1.0))
end
