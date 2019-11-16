using LavoisierCore
using Test

model = GERG2008(:C2,:C3)
x0 = [0.5,0.5]
T0 = 300
@testset "initial" begin
    @test length(model) == 2
end