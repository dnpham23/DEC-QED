
using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include("operators.jl");

# dimensions of the problem
###s = 1e-1;
W = 160.0;
Tm = 240.0;
dx = 0.02;
dt = 0.016;
q = 1.0;
eta = 0.002;
kappa = 0.006;
alpha = 0.0; # resistive loss
beta = 0.0; # normalized biased current density
sigon = 10.0;
sigoff = 10.0;
t_on = 3*sigon;
t_off= Tm-17*sigoff;
A = 1.4;
Omega = 0.6;

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
##xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
##xi0 = (xlist .- x0)/sqrt(1-u^2);

# rhs 
###rhs = zeros(nv_x);
rhs = -beta*ones(nv_x);

# jj current density
qq = q*ones(nv_x);


# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);

# initial conditions
###phit[1,:] = -dt*4*u*exp.(xi0)./(sqrt(1-u^2)*(1 .+ exp.(2*xi0)));
###phi[1,:]  = 4*atan.(exp.(xi0));
phi[1,:] .= asin(beta);
phit[1,:] .= 0.0;

# the on-off time indices
t_onend = Int(round(t_on/dt)) + 1
t_offstart = Int(round(t_off/dt)) + 1



for t = 1:nv_t-1
    
    phi_current  = phi[t,:]
    phix_current = phi[t,2:end] .- phi[t,1:end-1]
    phit_past    = phit[t,:]
    phix[t,:]    = phix_current # update the current phix 

    # compute the laplacian
    phit_current = sinegordon_dec_3(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, qq, alpha, "hard-wall");

    if t <= t_onend
        At = A*exp(-(t*dt-t_on)^2/(2*sigon^2))
    elseif t_onend < t < t_offstart
        At = A
    else
        At = A*exp(-(t*dt-t_off)^2/(2*sigoff^2))
    end

    phi[t+1,:]      = phit_current + phi_current
    phi[t+1,1]      = phi[t+1,2] - At*sin(Omega*t*dt)*dx
    phi[t+1,nv_x]   = phi[t+1,nv_x-1] # phix=0 at x=W

    phit_current[1] = phi[t+1,1] - phi[t,1]
    phit_current[nv_x] = phi[t+1,nv_x] - phi[t,nv_x]
    phit[t+1,:] = phit_current
end
phix[nv_t,:] = phi[nv_t,2:end] .- phi[nv_t,1:end-1]



using Plots
gr(size = (900, 500), legend = false, colorbar = true);

color_phixmax = findmax(phix)[1]/2;
color_phixmin = findmin(phix)[1]/2;
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
savefig(myplot,"phix.png")

color_phitmax = findmax(phit)[1]/3;
color_phitmin = findmin(phit)[1]/3;
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
savefig(myplot,"phit.png")

phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
savefig(myplot,"phi.png")




