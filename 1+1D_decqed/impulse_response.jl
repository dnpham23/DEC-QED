
using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include(joinpath(@__DIR__, "operators.jl"));
include(joinpath(@__DIR__, "simulation.jl"));


using CSV
using DataFrames

phi_ss = CSV.read("phiG1p5K0p2.csv", DataFrame);
phi_ss = collect(values(phi_ss[:,1]));
##phi_ss = phi_ss[1:end-1];
			
Nv_total_half = Int((length(phi_ss)-1)/2); 

# dimensions of the problem
dx = 0.01;
W = 20.0-2*dx; 
g = 1.5;
kappa = 0.2;
Tm = 20.0;
dt = 0.008;
nv_x_half = Int(round(W/dx/2));
nv_x = 2*nv_x_half+1;
#phi_ss = phi_ss[Nv_total_half+1-nv_x_half:Nv_total_half+1+nv_x_half];

V1  = g^2 .+ 2*pi*kappa*cos.(2*sqrt(pi).*phi_ss);

# dimensions of the problem
biasloc_start = 0.0;
biasloc_end   = W/2;
to = 40.0;
halfbandw = 16.0;
##halfbandw = 40.0;
nv_halfband = Int(round(halfbandw/dx));
W_zoom = 5.0;
nv_x_zoom = Int(round(W_zoom/dx)) + 1

# global nodal params
nv_x = Int(round(W/dx)) + 1
nv_t = Int(round(Tm/dt)) + 1
Nv   = nv_x*nv_t
xlist = collect(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx);

# rhs 
#biasloc_indstart = Int(round((biasloc_start-(-W/2))/dx)) + 1; 
#biasloc_indend = Int(round((biasloc_end-(-W/2))/dx)) + 1; 
rhs = zeros(nv_x);
rhs[nv_x_half+1] = 1.0/(dt*dx); # delta function at the center

# initial conditions
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);

operator! = (phi_current, phix_current, phit_past) ->
    linearized_fluctuation(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, V1, "radiative")
evolve!(phi, phit, phix, operator!)

phit_fin = phit[end,:];
phi_fin = phi[end,:];

##########################
# plotting

using Plots
gr(size = (1200, 900), legend = false, colorbar = true);
##gr(size = (500, 900), legend = false, colorbar = true);


color_phixmax = findmax(phix)[1];
color_phixmin = findmin(phix)[1];
myplot = heatmap(-(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, dt:dt:nv_t*dt, phix, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
savefig(myplot,"phix.png")
savefig(myplot,"phix.svg")

color_phitmax = findmax(phit)[1];
color_phitmin = findmin(phit)[1];
myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phit, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
savefig(myplot,"phit.png")
savefig(myplot,"phit.svg")

phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phi, xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
savefig(myplot,"phi.png")
savefig(myplot,"phi.svg")


# zooming into the soliton location
##gr(size = (400, 650), legend = false, colorbar = true);

#myplot = heatmap(-(nv_x_zoom-1)/2*dx+dx/2:dx:(nv_x_zoom-1)/2*dx-dx/2, dt:dt:nv_t*dt, phix[:,Int((nv_x-1)/2-(nv_x_zoom-1)/2+1) : Int((nv_x-1)/2+(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix_zoomed.png")

#myplot = heatmap(-(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, dt:dt:nv_t*dt, phit[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit_zoomed.png")

#myplot = heatmap(-(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, dt:dt:nv_t*dt, phi[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)], xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi_zoomed.png")


########################## 
###@save "W400_Tm150_u1_3_alpha0_beta0.jld2" potential kinetic jjkinetic energytotal phi W Tm dx dt q eta kappa u alpha beta x0 nv_x nv_t Nv # only saves phi for high-res calculations
###@save "energies_vav_W8000_Tm4000_u0_55_g0_4_t0=12_5.jld2" potential kinetic jjkinetic energyfield energytotal potential_bw kinetic_bw jjkinetic_bw energyfield_bw energytotal_bw W Tm dx dt q eta kappa u alpha beta x0 g nv_x nv_t Nv biasloc_start biasloc_end W_zoom nv_x_zoom halfbandw
###@save "endpoint_vav_W8000_Tm4000_u0_3_g0_325_t0=20.jld2" potential kinetic jjkinetic energyfield energytotal potential_bw kinetic_bw jjkinetic_bw energyfield_bw energytotal_bw W Tm dx dt q eta kappa u alpha beta x0 g nv_x nv_t Nv biasloc_start biasloc_end W_zoom nv_x_zoom halfbandw phit_fin phi_fin 
###@save "endpoint_single_W12000_Tm6000_u0_55_g0_28_x0=0.jld2" potential kinetic jjkinetic energyfield energytotal potential_bw kinetic_bw jjkinetic_bw energyfield_bw energytotal_bw W Tm dx dt q eta kappa u alpha beta x0 g nv_x nv_t Nv biasloc_start biasloc_end W_zoom nv_x_zoom phit_fin phi_fin 
###@save "W10_Tm10_kappa1_g1_5_dx0_01_dt0_008.jld2" W Tm dx dt kappa g nv_x nv_t Nv phit_fin phi_fin phi



##########################
### rotated 90 degrees ###
###gr(size = (5500, 420), legend = false, colorbar = true);
###gr(size = (3000, 450), legend = false, colorbar = true);
###gr(size = (2000, 500), legend = false, colorbar = true);
#gr(size = (1000, 600), legend = false, colorbar = true);


#color_phixmax = findmax(phix)[1];
#color_phixmin = findmin(phix)[1];
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix.png")

#color_phitmax = findmax(phit)[1];
#color_phitmin = findmin(phit)[1];
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit', ylabel =("x"), xlabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit.png")

#phimax = findmax(phi)[1];
#phimin = findmin(phi)[1];
#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi', ylabel =("x"), xlabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi.png")



## zooming into the soliton location
#gr(size = (2000, 400), legend = false, colorbar = true);

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx+dx/2:dx:(nv_x_zoom-1)/2*dx-dx/2, phix[:,Int((nv_x-1)/2-(nv_x_zoom-1)/2+1) : Int((nv_x-1)/2+(nv_x_zoom-1)/2)]' , xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (color_phixmin, color_phixmax));
#savefig(myplot,"phix_zoomed.png")

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, phit[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)]', xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (color_phitmin, color_phitmax));
#savefig(myplot,"phit_zoomed.png")

#myplot = heatmap(dt:dt:nv_t*dt, -(nv_x_zoom-1)/2*dx:dx:(nv_x_zoom-1)/2*dx, phi[:,Int((nv_x-1)/2+1-(nv_x_zoom-1)/2) : Int((nv_x-1)/2+1 +(nv_x_zoom-1)/2)]', xlabel =("x"), ylabel = ("t"), title = "\$\\phi\$", clims = (phimin, phimax));
#savefig(myplot,"phi_zoomed.png")



