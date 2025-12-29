function kleingordon_massive_dec(dx::Float64, dt::Float64, nv_x::Int64, phix_current::Array{Float64,1}, phit_past::Array{Float64,1},
                                 phi_current::Array{Float64,1}, rhs::Array{Float64,1}, g::Float64, m::Float64)

    phit_current = zeros(nv_x)

    da = dx*dt; # uniform grid
    tempx = dt/(dx*da)
    tempt = dx/(dt*da)
    for v = 2:nv_x-1 # exclude the boundary points
        phit_current[v] = ((phit_past[v]*tempt - phix_current[v-1]*tempx + phix_current[v]*tempx) - (rhs[v] + phi_current[v]*g^2 - pi*g*m*sin(2*sqrt(pi)*phi_current[v])) )/tempt
    end

    return phit_current
end


function sinegordon_dec(dx::Float64, dt::Float64, nv_x::Int64, phix_current::Array{Float64,1}, phit_past::Array{Float64,1},
                          phi_current::Array{Float64,1}, rhs::Array{Float64,1}, q::Float64, alpha::Float64)

    phit_current = zeros(nv_x)

    da = dx*dt; # uniform grid
    tempx = dt/(dx*da)
    tempr1 = 1+alpha*dt/2
    tempr2 = 1-alpha*dt/2
    for v = 2:nv_x-1 # exclude the boundary points
        phit_current[v] = tempr2/tempr1*phit_past[v] + dt^2/tempr1*( (phix_current[v]-phix_current[v-1])/dx^2 - rhs[v] - q*sin(phi_current[v]) )
    end

    return phit_current
end

function sinegordon_dec_2(dx::Float64, dt::Float64, nv_x::Int64, phix_current::Array{Float64,1}, phit_past::Array{Float64,1},
                          phi_current::Array{Float64,1}, rhs::Array{Float64,1}, q::Float64, alpha::Float64, bctype::String)

    phit_current = zeros(nv_x)

    da = dx*dt; # uniform grid
    tempx = dt/(dx*da)
    tempr1 = 1+alpha*dt/2
    tempr2 = 1-alpha*dt/2
    for v = 2:nv_x-1 # exclude the boundary points
        phit_current[v] = tempr2/tempr1*phit_past[v] + dt^2/tempr1*( (phix_current[v]-phix_current[v-1])/dx^2 - rhs[v] - q*sin(phi_current[v]) )
    end 
    
    if bctype == "periodic"
       phit_current[1] = tempr2/tempr1*phit_past[1] + dt^2/tempr1*( (phix_current[1]-phix_current[nv_x-1])/dx^2 - rhs[1] - q*sin(phi_current[1]) )
       phit_current[nv_x] = phit_current[1]
    end

    return phit_current
end

function sinegordon_dec_3(dx::Float64, dt::Float64, nv_x::Int64, phix_current::Array{Float64,1}, phit_past::Array{Float64,1},
                          phi_current::Array{Float64,1}, rhs::Array{Float64,1}, qq::Array{Float64,1}, alpha::Float64, bctype::String)

    phit_current = zeros(nv_x)

    da = dx*dt; # uniform grid
    tempx = dt/(dx*da)
    tempr1 = 1+alpha*dt/2
    tempr2 = 1-alpha*dt/2
    for v = 2:nv_x-1 # exclude the boundary points
        phit_current[v] = tempr2/tempr1*phit_past[v] + dt^2/tempr1*( (phix_current[v]-phix_current[v-1])/dx^2 - rhs[v] - qq[v]*sin(phi_current[v]) )
    end 
    
    if bctype == "periodic"
       phit_current[1] = tempr2/tempr1*phit_past[1] + dt^2/tempr1*( (phix_current[1]-phix_current[nv_x-1])/dx^2 - rhs[1] - qq[1]*sin(phi_current[1]) )
       phit_current[nv_x] = phit_current[1]
    end

    return phit_current
end