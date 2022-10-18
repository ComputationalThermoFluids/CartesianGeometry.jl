using CartesianGeometry
using Test

@testset "CartesianGeometry.jl" begin
    range = 1:8

    x = collocated(identity, range, range)
    @test first(x) == zero(eltype(x))
    @test last(x) == one(eltype(x))

    x = staggered(identity, range, range)
    @test 2length(x) * first(x) == one(eltype(x))
end
