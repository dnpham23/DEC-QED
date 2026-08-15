"""Finite-difference DEC evolution operators for the 1+1D examples."""

function kleingordon_massive_dec(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                                 phi_current::AbstractVector, rhs::AbstractVector, g::Real, m::Real)
    phit_current = zeros(nv_x)
    for v = 2:nv_x-1
        phit_current[v] = phit_past[v] + dt^2 * ((phix_current[v] - phix_current[v - 1]) / dx^2 - rhs[v] - phi_current[v] * g^2 + pi * g * m * sin(2 * sqrt(pi) * phi_current[v]))
    end
    return phit_current
end

function sinegordon_dec(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                        phi_current::AbstractVector, rhs::AbstractVector, q::Real, alpha::Real)
    phit_current = zeros(nv_x)
    temporal_scale = (1 - alpha * dt / 2) / (1 + alpha * dt / 2)
    spatial_scale = dt^2 / (1 + alpha * dt / 2)
    for v = 2:nv_x-1
        phit_current[v] = temporal_scale * phit_past[v] + spatial_scale * ((phix_current[v] - phix_current[v - 1]) / dx^2 - rhs[v] - q * sin(phi_current[v]))
    end
    return phit_current
end

function sinegordon_dec_2(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                          phi_current::AbstractVector, rhs::AbstractVector, q::Real, alpha::Real, bctype::AbstractString)
    phit_current = sinegordon_dec(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, q, alpha)
    if bctype == "periodic"
        temporal_scale = (1 - alpha * dt / 2) / (1 + alpha * dt / 2)
        spatial_scale = dt^2 / (1 + alpha * dt / 2)
        phit_current[1] = temporal_scale * phit_past[1] + spatial_scale * ((phix_current[1] - phix_current[end]) / dx^2 - rhs[1] - q * sin(phi_current[1]))
        phit_current[end] = phit_current[1]
    end
    return phit_current
end

function sinegordon_dec_3(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                          phi_current::AbstractVector, rhs::AbstractVector, qq::AbstractVector, alpha::Real, bctype::AbstractString)
    phit_current = zeros(nv_x)
    temporal_scale = (1 - alpha * dt / 2) / (1 + alpha * dt / 2)
    spatial_scale = dt^2 / (1 + alpha * dt / 2)
    for v = 2:nv_x-1
        phit_current[v] = temporal_scale * phit_past[v] + spatial_scale * ((phix_current[v] - phix_current[v - 1]) / dx^2 - rhs[v] - qq[v] * sin(phi_current[v]))
    end
    if bctype == "periodic"
        phit_current[1] = temporal_scale * phit_past[1] + spatial_scale * ((phix_current[1] - phix_current[end]) / dx^2 - rhs[1] - qq[1] * sin(phi_current[1]))
        phit_current[end] = phit_current[1]
    end
    return phit_current
end

function laplacian_dec_spacetime(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                                 phi_current::AbstractVector, rhs::AbstractVector, g::Real)
    phit_current = zeros(nv_x)
    for v = 2:nv_x-1
        phit_current[v] = phit_past[v] + dt^2 * ((phix_current[v] - phix_current[v - 1]) / dx^2 - rhs[v] - phi_current[v] * g^2)
    end
    return phit_current
end

function sinekleingordon_dec(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                             phi_current::AbstractVector, rhs::AbstractVector, qq::AbstractVector, alpha::Real, bctype::AbstractString, g::Real)
    phit_current = sinegordon_dec_3(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, qq, alpha, bctype)
    mass_scale = dt^2 / (1 + alpha * dt / 2)
    for v = 2:nv_x-1
        phit_current[v] -= mass_scale * g^2 * phi_current[v]
    end
    if bctype == "periodic"
        phit_current[1] -= mass_scale * g^2 * phi_current[1]
        phit_current[end] = phit_current[1]
    end
    return phit_current
end

function linearized_fluctuation(dx::Real, dt::Real, nv_x::Integer, phix_current::AbstractVector, phit_past::AbstractVector,
                                phi_current::AbstractVector, rhs::AbstractVector, V1::AbstractVector, bctype::AbstractString)
    phit_current = zeros(nv_x)
    for v = 2:nv_x-1
        phit_current[v] = phit_past[v] + dt^2 * ((phix_current[v] - phix_current[v - 1]) / dx^2 + rhs[v] - V1[v] * phi_current[v])
    end
    if bctype == "periodic"
        phit_current[1] = phit_past[1] + dt^2 * ((phix_current[1] - phix_current[end]) / dx^2 + rhs[1] - V1[1] * phi_current[1])
        phit_current[end] = phit_current[1]
    elseif bctype == "radiative"
        phit_current[1] = phix_current[1] * dt / dx
        phit_current[end] = -phix_current[end] * dt / dx
    end
    return phit_current
end
