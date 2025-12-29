
using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include("operators.jl");

# dimensions of the problem
###s = 1e-1;
W = 100.0;
Tm = 1000.0;
dx = 0.02;
dt = 0.016;
q = 1.0;
eta = 0.002;
kappa = 0.006;
u = 0.55;
alpha = 0.003; # resistive loss
beta = 0.001; # normalized biased current density

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
xi0 = xlist/sqrt(1-u^2);

# rhs 
###rhs = zeros(nv_x);
rhs = -beta*ones(nv_x);


# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);
# initial conditions
phit[1,:] = -dt*4*u*exp.(xi0)./(sqrt(1-u^2)*(1 .+ exp.(2*xi0)));
phi[1,:]  = 4*atan.(exp.(xi0));

for t = 1:nv_t-1
    phi_current  = phi[t,:]
    phix_current = phi[t,2:end] .- phi[t,1:end-1]
    phit_past    = phit[t,:]
    phix[t,:]    = phix_current # update the current phix 

    # compute the laplacian
    phit_current = sinegordon_dec_2(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, q, alpha, "periodic");
    phi[t+1,:]   = phit_current + phi_current 
    phit[t+1,:]  = phit_current
end
phix[nv_t,:] = phi[nv_t,2:end] .- phi[nv_t,1:end-1]



using Plots
gr(size = (3000, 450), legend = false, colorbar = true);

color_phix = findmax(phix)[1]/2;
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (-color_phix, color_phix));
savefig(myplot,"phix.png")

color_phit = findmax(phit)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (-color_phit, color_phit));
savefig(myplot,"phit.png")

phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
savefig(myplot,"phi.png")




