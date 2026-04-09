
using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include("operators.jl");

# dimensions of the problem
###s = 1e-1;
W = 300.0;
Tm = 24800.0;
dx = 0.05;
dt = 0.04;
###q = -1.0; # negative sine term
q = 1.0; 
eta = 0.0;
kappa = 0.0;
###u = 1.3; # u >1 
u = 0.55;
alpha = 0.0; # resistive loss
g = 0.3;
beta = 2*pi*g^2; # to ensure zero total charges, beta depends on g
###beta = 0.0;
x0 = -10.0; # starting point of the fluxon 
biasloc_start = 0.0;
biasloc_end   = W/2;
to = 40.0;
halfbandw = 16.0;
##halfbandw = 40.0;
nv_halfband = Int(round(halfbandw/dx));

W_zoom = 30.0;
##W_zoom = 80.0;
nv_x_zoom = Int(round(W_zoom/dx)) + 1


# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);
gamma = 1/sqrt(1-u^2)
###xi0 = -(xlist .- x0)/sqrt(u^2-1); # opposite sign to the original SG
xi0 = (xlist .- x0)/sqrt(1-u^2);


# rhs 
###rhs = -beta*ones(nv_x);
biasloc_indstart = Int(round((biasloc_start-(-W/2))/dx)) + 1; 
biasloc_indend = Int(round((biasloc_end-(-W/2))/dx)) + 1; 
rhs = zeros(nv_x);
rhs[biasloc_indstart:biasloc_indend] .= -beta;


# jj current density
qq = q*ones(nv_x); 

# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);

# initial conditions
###phit[1,:] = dt*4*u*exp.(xi0)./(sqrt(u^2-1)*(1 .+ exp.(2*xi0))); # opposite sign to the original SG
phit[1,:] = -dt*4*u*exp.(xi0)./(sqrt(1-u^2)*(1 .+ exp.(2*xi0)));

phi[1,:]  = 4*atan.(exp.(xi0)); # OG SG fluxon solution

# fluxon-antifluxon pair
###phit[1,:] = -dt*4*gamma./( 1 .+ ( sinh.(-u*to*gamma)./(u.*cosh.((xlist.-x0).*gamma)) ).^2 ).*cosh.(-u*to*gamma)./cosh.((xlist.-x0).*gamma);
###phi[1,:] = -4*atan.( sinh.(-u*to*gamma)./(u.*cosh.((xlist.-x0).*gamma)) );

phit[1,1:Int((nv_x-1)/2)-5000] .= 0.0;
phit[1,Int((nv_x-1)/2)+5000:end] .= 0.0;



for t = 1:nv_t-1
    phi_current  = phi[t,:]
    phix_current = phi[t,2:end] .- phi[t,1:end-1]
    phit_past    = phit[t,:]
    phix[t,:]    = phix_current # update the current phix 

    # compute the laplacian
    phit_current = sinekleingordon_dec(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, qq, alpha, "hard-wall", g);
    phi[t+1,:]      = phit_current + phi_current 
    phi[t+1,1]      = phi[t+1,2] - (eta+kappa)*dx
    phi[t+1,nv_x]   = phi[t+1,nv_x-1] + (eta-kappa)*dx
    phit_current[1] = phi[t+1,1] - phi[t,1]
    phit_current[nv_x] = phi[t+1,nv_x] - phi[t,nv_x]
    phit[t+1,:] = phit_current
end
phix[nv_t,:] = phi[nv_t,2:end] .- phi[nv_t,1:end-1]

phit_fin = phit[end,:];
phi_fin = phi[end,:];


###########################
# compute the energies

potential = 0.5/dx*sum(phix.^2, dims=2);
kinetic   = 0.5*dx/dt^2*sum(phit.^2, dims=2);
jjkinetic = dx*sum(1 .- q.*cos.(phi), dims=2);

extfield = zeros(nv_x)';
extfield[biasloc_indstart:biasloc_indend] .= -beta/g;
efield = g.*phi .+ extfield;
energyfield = 0.5*dx*sum(efield.^2, dims=2);

energytotal = potential[1:end-1] + jjkinetic[1:end-1] + 0.5*(kinetic[1:end-1]+kinetic[2:end]) + energyfield[1:end-1] 


# energy within a finite region around the origin

potential_bw = 0.5/dx*sum(phix[:,Int((nv_x-1)/2)+1-nv_halfband:Int((nv_x-1)/2)+nv_halfband].^2, dims=2);
kinetic_bw   = 0.5*dx/dt^2*sum(phit[:,Int((nv_x-1)/2)+1-nv_halfband:Int((nv_x-1)/2)+1+nv_halfband].^2, dims=2);
jjkinetic_bw = dx*sum(1 .- q.*cos.(phi[:,Int((nv_x-1)/2)+1-nv_halfband:Int((nv_x-1)/2)+1+nv_halfband]), dims=2);

energyfield_bw = 0.5*dx*sum(efield[:,Int((nv_x-1)/2)+1-nv_halfband:Int((nv_x-1)/2)+1+nv_halfband].^2, dims=2);

energytotal_bw = potential_bw[1:end-1] + jjkinetic_bw[1:end-1] + 0.5*(kinetic_bw[1:end-1]+kinetic_bw[2:end]) + energyfield_bw[1:end-1] 



##########################
# plotting

using Plots
#gr(size = (1500, 650), legend = false, colorbar = true);
##gr(size = (500, 900), legend = false, colorbar = true);


#color_phixmax = findmax(phix)[1]/2;
#color_phixmin = findmin(phix)[1]/2;
#myplot = heatmap(-(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, dt:dt:nv_t*dt, phix, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix.png")
###savefig(myplot,"phix.svg")

#color_phitmax = findmax(phit)[1]/3;
#color_phitmin = findmin(phit)[1]/3;
#myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phit, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit.png")
###savefig(myplot,"phit.svg")

#phimax = findmax(phi)[1];
#phimin = findmin(phi)[1];
#myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phi, xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi.png")
###savefig(myplot,"phi.svg")


# zooming into the soliton location
##gr(size = (400, 650), legend = false, colorbar = true);

#myplot = heatmap(-(nv_x_zoom-1)/2*dx+dx/2:dx:(nv_x_zoom-1)/2*dx-dx/2, dt:dt:nv_t*dt, phix[:,Int((nv_x-1)/2-(nv_x_zoom-1)/2+1) : Int((nv_x-1)/2+(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix_zoomed.png")

#myplot = heatmap(-(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, dt:dt:nv_t*dt, phit[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit_zoomed.png")

#myplot = heatmap(-(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, dt:dt:nv_t*dt, phi[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi_zoomed.png")


########################## 
@save "single_W300_Tm24800_u0_55_g0_3_x0=-10_dx0_05_dt0_04.jld2" potential kinetic jjkinetic energyfield energytotal potential_bw kinetic_bw jjkinetic_bw energyfield_bw energytotal_bw W Tm dx dt q eta kappa u alpha beta x0 g nv_x nv_t Nv biasloc_start biasloc_end W_zoom nv_x_zoom phit_fin phi_fin


##########################
### rotated 90 degrees ###
gr(size = (1000, 600), legend = false, colorbar = true);


#color_phixmax = findmax(phix)[1]/2;
#color_phixmin = findmin(phix)[1]/2;
####color_phix = findmax(phix)[1]/8;
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix.png")

###color_phitmax = findmax(phit)[1]*3;
#color_phitmax = findmax(phit)[1]/4;
#color_phitmin = findmin(phit)[1]/4;
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit.png")

#phimax = findmax(phi)[1];
#phimin = findmin(phi)[1];
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi.png")



## zooming into the soliton location
gr(size = (2000, 400), legend = false, colorbar = true);

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx+dx/2:dx:(nv_x_zoom-1)/2*dx-dx/2, phix[:,Int((nv_x-1)/2-(nv_x_zoom-1)/2+1) : Int((nv_x-1)/2+(nv_x_zoom-1)/2)]' , xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix_zoomed.png")

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, phit[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)]', xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit_zoomed.png")

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, phi[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)]', xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi_zoomed.png")



