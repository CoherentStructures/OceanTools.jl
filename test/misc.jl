using OceanTools, Test

@testset "misc" begin
    lat = rand()
    @test OceanTools.earthDiffTensor([0., lat]) ≈ [1/cosd(lat)^2 0; 0 1]
end