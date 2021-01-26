"""
    earthDiffTensor(x)

Compute the isotropic diffusion tensor at location `x`, given in lon-lat
coordinates (in degrees).
"""
earthDiffTensor(x) = Tensors.SymmetricTensor{2,2}((1/cosd(x[2])^2, 0, 1))
