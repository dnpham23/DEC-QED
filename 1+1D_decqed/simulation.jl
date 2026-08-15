"""
Shared utilities for the 1+1D DEC-QED examples.

The functions in this file deliberately have no package dependencies, so examples can
be run from any working directory after including this file.
"""

"""Build a centered, uniformly spaced spatial grid."""
function centered_grid(width::Real, dx::Real)
    width > 0 || throw(ArgumentError("width must be positive"))
    dx > 0 || throw(ArgumentError("dx must be positive"))

    nv_x = round(Int, width / dx) + 1
    x = collect(range(-(nv_x - 1) * dx / 2, step = dx, length = nv_x))
    return nv_x, x
end

"""Return the number of time nodes for a uniform time grid."""
function time_nodes(duration::Real, dt::Real)
    duration > 0 || throw(ArgumentError("duration must be positive"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    return round(Int, duration / dt) + 1
end

"""Apply a fixed spatial-gradient boundary condition to a field increment."""
function fixed_gradient_boundary!(next_phi, current_phi, increment, _, left_gradient, right_gradient, dx)
    next_phi[1] = next_phi[2] - left_gradient * dx
    next_phi[end] = next_phi[end - 1] + right_gradient * dx
    increment[1] = next_phi[1] - current_phi[1]
    increment[end] = next_phi[end] - current_phi[end]
    return nothing
end

"""
    evolve!(phi, phit, phix, operator!; boundary! = nothing)

Advance preallocated field histories in place. `operator!` receives `(phi, phix,
phit)` for the current time step and returns the field increment for the next step.
`boundary!`, when supplied, receives `(next_phi, current_phi, increment, step)` and
must keep `increment` consistent with any boundary values it changes.
"""
function evolve!(phi::AbstractMatrix, phit::AbstractMatrix, phix::AbstractMatrix, operator!; boundary! = nothing)
    nv_t, nv_x = size(phi)
    nv_x >= 3 || throw(ArgumentError("at least three spatial nodes are required"))
    size(phit) == (nv_t, nv_x) || throw(DimensionMismatch("phit must match phi"))
    size(phix) == (nv_t, nv_x - 1) || throw(DimensionMismatch("phix must have one fewer spatial column than phi"))

    @views for step in 1:nv_t-1
        current_phi = phi[step, :]
        current_phit = phit[step, :]
        current_phix = phix[step, :]

        @inbounds for x in eachindex(current_phix)
            current_phix[x] = current_phi[x + 1] - current_phi[x]
        end

        increment = operator!(current_phi, current_phix, current_phit)
        length(increment) == nv_x || throw(DimensionMismatch("operator returned an increment with the wrong length"))

        next_phi = phi[step + 1, :]
        next_phit = phit[step + 1, :]
        @. next_phi = current_phi + increment
        boundary! !== nothing && boundary!(next_phi, current_phi, increment, step)
        next_phit .= increment
    end

    @views begin
        final_phi = phi[end, :]
        final_phix = phix[end, :]
        @inbounds for x in eachindex(final_phix)
            final_phix[x] = final_phi[x + 1] - final_phi[x]
        end
    end
    return nothing
end