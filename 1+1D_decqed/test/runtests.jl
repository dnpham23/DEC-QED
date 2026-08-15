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

    @test_throws ArgumentError centered_grid(0.0, 1.0)
    @test_throws ArgumentError time_nodes(1.0, 0.0)
    @test_throws ArgumentError evolve!(zeros(2, 2), zeros(2, 2), zeros(2, 1), operator!)
    @test_throws DimensionMismatch evolve!(zeros(2, 3), zeros(2, 2), zeros(2, 2), operator!)
    @test_throws DimensionMismatch evolve!(zeros(2, 3), zeros(2, 3), zeros(2, 2), (_, _, _) -> zeros(2))
end

@testset "Klein-Gordon operators" begin
    phi = [0.5, 1.0, 1.5, 2.0]
    phix = [0.25, 0.75, 1.25]
    phit = [0.1, 0.2, 0.3, 0.4]
    rhs = [0.0, 0.2, -0.1, 0.0]

    massless = laplacian_dec_spacetime(0.5, 0.25, 4, phix, phit, phi, rhs, 2.0)
    @test massless == [0.0, 0.2 + 0.25^2 * ((0.75 - 0.25) / 0.5^2 - 0.2 - 2.0^2),
                        0.3 + 0.25^2 * ((1.25 - 0.75) / 0.5^2 + 0.1 - 2.0^2 * 1.5), 0.0]

    massive = kleingordon_massive_dec(0.5, 0.25, 4, phix, phit, phi, rhs, 2.0, 0.3)
    @test massive[1] == massive[end] == 0.0
    @test massive[2] ≈ phit[2] + 0.25^2 * ((phix[2] - phix[1]) / 0.5^2 - rhs[2] - phi[2] * 2.0^2 + pi * 2.0 * 0.3 * sin(2 * sqrt(pi) * phi[2]))
end

@testset "Sine-Gordon operators" begin
    phi = [1.0, 0.2, -0.4, 0.5]
    phix = [0.3, -0.1, 0.2]
    phit = [0.1, 0.2, 0.3, 0.4]
    rhs = [0.4, 0.5, 0.6, 0.7]
    qq = [1.0, 1.5, 2.0, 2.5]

    base = sinegordon_dec(0.5, 0.25, 4, phix, phit, phi, rhs, 1.2, 0.4)
    hard_wall = sinegordon_dec_2(0.5, 0.25, 4, phix, phit, phi, rhs, 1.2, 0.4, "hard-wall")
    @test hard_wall == base

    periodic = sinegordon_dec_2(0.5, 0.25, 4, phix, phit, phi, rhs, 1.2, 0.4, "periodic")
    @test periodic[end] == periodic[1]
    @test periodic[2:3] == base[2:3]

    variable_q = sinegordon_dec_3(0.5, 0.25, 4, phix, phit, phi, rhs, qq, 0.4, "periodic")
    @test variable_q[end] == variable_q[1]

    sine_klein = sinekleingordon_dec(0.5, 0.25, 4, phix, phit, phi, rhs, qq, 0.4, "periodic", 2.0)
    mass_scale = 0.25^2 / (1 + 0.4 * 0.25 / 2)
    @test sine_klein[1] ≈ variable_q[1] - mass_scale * 2.0^2 * phi[1]
    @test sine_klein[2] ≈ variable_q[2] - mass_scale * 2.0^2 * phi[2]
    @test sine_klein[end] == sine_klein[1]
end

@testset "Linearized fluctuation boundaries" begin
    phi = [1.0, 0.0, 0.0, 0.0]
    phix = [0.3, 0.0, -0.2]
    phit = zeros(4)
    rhs = zeros(4)
    potential = zeros(4)

    periodic = linearized_fluctuation(0.5, 0.25, 4, phix, phit, phi, rhs, potential, "periodic")
    @test periodic[end] == periodic[1]

    radiative = linearized_fluctuation(0.5, 0.25, 4, phix, phit, phi, rhs, potential, "radiative")
    @test radiative[1] == 0.3 * 0.25 / 0.5
    @test radiative[end] == -(-0.2) * 0.25 / 0.5
end