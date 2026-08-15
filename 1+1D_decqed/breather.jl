
using JLD2;
using SparseArrays
using LinearAlgebra
include(joinpath(@__DIR__, "operators.jl"));
include(joinpath(@__DIR__, "simulation.jl"));

# dimensions of the problem
W = 100.0;
Tm = 1500.0;
dx = 0.02;
dt = 0.016;
q = 1.0;
eta = 0.002;
kappa = 0.006;
u = 0.55;
alpha = 0.0; # resistive loss
beta = 0.0; # normalized biased current density
nu = 0.6;

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
gamma = 1/sqrt(1-u^2); 
xi0 = xlist*gamma;

# rhs 
###rhs = zeros(nv_x);
rhs = -beta*ones(nv_x);


# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);
# initial conditions
phit[1,:] = dt*4*tan(nu)./(1 .+ (tan(nu)*sin.(-cos(nu)*u.*xi0) ./ cosh.(sin(nu).*xi0) ).^2).*( (cos(nu)*gamma*cos.(-cos(nu)*u.*xi0) ) ./ cosh.(sin(nu).*xi0) + (sin(nu)*gamma*u*sin.(-cos(nu)*u*xi0).*sinh.(sin(nu).*xi0))./(cosh.(sin(nu).*xi0)).^2 );
phi[1,:]  = 4*atan.(tan(nu)*sin.(-cos(nu)*u.*xi0) ./ cosh.(sin(nu).*xi0));

operator! = (phi_current, phix_current, phit_past) ->
    sinegordon_dec_2(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, q, alpha, "hard-wall")
boundary! = (next_phi, current_phi, increment, step) ->
    fixed_gradient_boundary!(next_phi, current_phi, increment, step, eta + kappa, eta - kappa, dx)
evolve!(phi, phit, phix, operator!; boundary! = boundary!)



using Plots
gr(size = (3000, 420), legend = false, colorbar = true);


color_phix = findmax(phix)[1]/2;
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (-color_phix, color_phix));
savefig(myplot,"phix.png")

color_phit = findmax(phit)[1]/2;
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (-color_phit, color_phit));
savefig(myplot,"phit.png")

phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
savefig(myplot,"phi.png")




