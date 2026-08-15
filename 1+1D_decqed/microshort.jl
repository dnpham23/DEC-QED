
using JLD2;
using SparseArrays
using LinearAlgebra
include(joinpath(@__DIR__, "operators.jl"));
include(joinpath(@__DIR__, "simulation.jl"));

# dimensions of the problem
W = 40.0;
Tm = 1500.0;
dx = 0.02;
dt = 0.016;
q = 1.0;
eta = 0.002;
kappa = 0.006;
alpha = 0.005; # resistive loss
beta = 0.002; # normalized biased current density
xs = -10.0; # location of the short
mu = 0.5; # JJ current of the short 
x0 = 10.0; # starting point of the fluxon 
u = -1/sqrt(1+(4*alpha/(pi*beta))^2);

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
xi0 = (xlist .- x0)/sqrt(1-u^2);

# rhs 
###rhs = zeros(nv_x);
rhs = -beta*ones(nv_x);

# jj current density
qq = q*ones(nv_x);
shortloc = Int(round((xs-(-W/2))/dx)) + 1; 
qq[shortloc] = q + mu/dx;


# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);
# initial conditions
phit[1,:] = -dt*4*u*exp.(xi0)./(sqrt(1-u^2)*(1 .+ exp.(2*xi0)));
phi[1,:]  = 4*atan.(exp.(xi0));

operator! = (phi_current, phix_current, phit_past) ->
    sinegordon_dec_3(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, qq, alpha, "hard-wall")
boundary! = (next_phi, current_phi, increment, step) ->
    fixed_gradient_boundary!(next_phi, current_phi, increment, step, eta + kappa, eta - kappa, dx)
evolve!(phi, phit, phix, operator!; boundary! = boundary!)



using Plots
gr(size = (3000, 450), legend = false, colorbar = true);

color_phix = findmax(phix)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (-color_phix, color_phix));
savefig(myplot,"phix.png")

color_phitmax = findmax(phit)[1];
color_phitmin = findmin(phit)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
savefig(myplot,"phit.png")

phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
savefig(myplot,"phi.png")




