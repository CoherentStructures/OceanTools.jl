using Coverage
import Pkg
Pkg.add("ForwardDiff")

cd(joinpath(@__DIR__, "..", "..")) do
    Codecov.submit(Codecov.process_folder())
end
