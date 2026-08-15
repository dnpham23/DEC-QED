using Test

include(joinpath(@__DIR__, "..", "operators.jl"))
include(joinpath(@__DIR__, "..", "simulation.jl"))

@testset "1+1D simulation utilities" begin
    nv_x, x = centered_grid(4.0, 1.0)
    @test nv_x == 5
    @test x == [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test time_nodes(1.0, 0.25) == 5

    phi = zeros(3, 5)
    phit = zeros(3, 5)
    phix = zeros(3, 4)
    operator! = (current_phi, current_phix, current_phit) -> ones(5)
    boundary! = (next_phi, current_phi, increment, step) ->
        fixed_gradient_boundary!(next_phi, current_phi, increment, step, 2.0, -1.0, 0.5)

    evolve!(phi, phit, phix, operator!; boundary! = boundary!)
    @test phi[2, :] == [0.0, 1.0, 1.0, 1.0, 0.5]
    @test phit[2, :] == [0.0, 1.0, 1.0, 1.0, 0.5]
    @test phix[end, :] == [1.0, 0.0, 0.0, -0.5]
end

@testset "periodic sine-Klein-Gordon boundary" begin
    phi = [1.0, 0.0, 0.0, 0.0]
    phix = [0.0, 0.0, 0.0]
    phit = zeros(4)
    rhs = zeros(4)
    qq = zeros(4)

    increment = sinekleingordon_dec(1.0, 1.0, 4, phix, phit, phi, rhs, qq, 0.0, "periodic", 2.0)
    @test increment[1] == -4.0
    @test increment[end] == increment[1]
end