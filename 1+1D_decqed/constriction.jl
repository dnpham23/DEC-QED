
using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include("operators.jl");

# dimensions of the problem
###s = 1e-1;
W = 200.0;
Tm = 1000.0;
dx = 0.02;
dt = 0.016;
q = 1.0;
eta = 0.002;
kappa = 0.006;
u = 0.85;
alpha = 0.0; # resistive loss
beta = 0.0; # normalized biased current density
xo = -50.0;
l = 20.0;
d = 30.0;
Jcs = 10.0;

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
xi0 = (xlist .- xo)/sqrt(1-u^2);

# rhs 
###rhs = zeros(nv_x);
rhs = -beta*ones(nv_x);

# the tapered regions
taperleft_startind = Int(round((-d-(-W/2))/dx)) + 1; 
taperleft_endind = Int(round((-l-(-W/2))/dx)) + 1; 
taperright_startind = Int(round((l-(-W/2))/dx)) + 1; 
taperright_endind = Int(round((d-(-W/2))/dx)) + 1; 


# jj current density
qq = q*ones(nv_x);
qq[taperleft_endind+1:taperright_startind-1] .= Jcs;
qq[taperright_startind:taperright_endind] = (d-l)./((1-1/Jcs).*(l:dx:d) .+ d*1/Jcs .- l);
qq[taperleft_startind:taperleft_endind] = reverse((d-l)./((1-1/Jcs).*(l:dx:d) .+ d*1/Jcs .- l)); 

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
    phit_current = sinegordon_dec_3(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, qq, alpha, "hard-wall");
    phi[t+1,:]      = phit_current + phi_current 
    phi[t+1,1]      = phi[t+1,2] - (eta+kappa)*dx
    phi[t+1,nv_x]   = phi[t+1,nv_x-1] + (eta-kappa)*dx
    phit_current[1] = phi[t+1,1] - phi[t,1]
    phit_current[nv_x] = phi[t+1,nv_x] - phi[t,nv_x]
    phit[t+1,:] = phit_current
end
phix[nv_t,:] = phi[nv_t,2:end] .- phi[nv_t,1:end-1]



using Plots
gr(size = (3000, 450), legend = false, colorbar = true);
###gr(size = (2000, 500), legend = false, colorbar = true);


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




