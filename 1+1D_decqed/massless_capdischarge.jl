#function main()

using JLD2;
#using IterativeSolvers;
using SparseArrays
using LinearAlgebra
include("operators.jl");

# dimensions of the problem
W = 6000.0;
g = 1.2;
L = 40.0;
Tm = 3000.0;
dx = 0.05;
dt = 0.04;
Q = 4.0;
capC = g*Q;


# global nodal params
nv_x = Int(W/dx) + 1
nv_t = Int(Tm/dt) + 1
Nv   = nv_x*nv_t

# find the locations of the capacitors 
mplane = Int((nv_x-1)/2)+1 - Int(L/(2*dx));
pplane = Int((nv_x-1)/2)+1 + Int(L/(2*dx));

# rhs 
rhs = zeros(nv_x);
rhs[mplane:pplane] .= capC;
#rhs[mplane+1:pplane-1] .= capC;
#rhs[mplane] = capC/2;
#rhs[pplane] = capC/2;

# time evolve
phit = zeros(nv_t, nv_x);
phix = zeros(nv_t, nv_x-1);
phi = zeros(nv_t, nv_x);
# initial conditions
phit[1,:] = zeros(nv_x);
phi[1,:]  = zeros(nv_x);
###phix[1,:] = zeros(nv_x-1); # redundant

for t = 1:nv_t-1
    phi_current  = phi[t,:]
    phix_current = phi[t,2:end] .- phi[t,1:end-1]
    phit_past    = phit[t,:]
    phix[t,:]    = phix_current # update the current phix 

    # compute the laplacian
    phit_current = laplacian_dec_spacetime(dx, dt, nv_x, phix_current, phit_past, phi_current, rhs, g);
    phit[t+1,:] = phit_current
    phi[t+1,:] = phit_current + phi_current 
end
phix[nv_t,:] = phi[nv_t,2:end] .- phi[nv_t,1:end-1];


currrent_charge_conserve_check = -phix[2:nv_t,:] + phix[1:nv_t-1,:] + phit[2:nv_t,2:nv_x] - phit[2:nv_t, 1:nv_x-1];

using Plots
#gr(size = (1500, 650), legend = false, colorbar = true);

#myplot = heatmap(-(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, dt:dt:(nv_t-1)*dt, currrent_charge_conserve_check, xlabel =("x"), ylabel = ("t") , title = "\$\\Phi(c)/\\phi_x \$");
#savefig(myplot,"chargeconserve_check.png")

#color_phix = findmax(phix)[1]/5;
#myplot = heatmap(-(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, dt:dt:nv_t*dt, phix, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_x\$", clims = (-color_phix, color_phix));
#savefig(myplot,"phix.png")

#color_phit = findmax(phit)[1]/5;
#myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phit, xlabel =("x"), ylabel = ("t"), title = "\$\\phi_t\$", clims = (-color_phit, color_phit));
#savefig(myplot,"phit.png")

#color_phi = findmax(phi)[1];
#myplot = heatmap(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, dt:dt:nv_t*dt, phi, xlabel =("x"), ylabel = ("t"), title = "\$\\psi\$", clims = (-color_phi, color_phi));
#savefig(myplot,"phi.png")



# calculating the  energies
energycharge = 0.5/dx*sum(phix.^2, dims=2);
energycurrent = 0.5*dx/dt^2*sum(phit.^2, dims=2);

extfield = zeros(nv_x)';
extfield[mplane:pplane] .= Q;
efield = g.*phi .+ extfield;

energyfield = 0.5*dx*sum(efield.^2, dims=2);


########################## 
###@save "SimulationData.jld2" phi W g L Tm dx dt Q nv_x nv_t Nv energycharge energycurrent energyfield # only saves phi for high-res calculations


using Plots
gr(size = (1500, 600), legend = false, colorbar = false)

startcut = 1;
endcut = nv_t-1;
x_loc = pplane + (pplane-mplane)*2
tempphi_vstime_list = phi[startcut:endcut,x_loc];
t = startcut*dt:dt:endcut*dt
myplot = plot(t,tempphi_vstime_list, title = "\$\\phi(t)\$ at \$x=100\$", xlabel="t", ylabel="\$\\phi\$");
savefig(myplot,"phi_vs_t_at_x=100_Tm3000.png");

temp_phit_vstime_list = phit[startcut:endcut,x_loc];
myplot = plot(t,temp_phit_vstime_list, title = "\$\\phi_t(t)\$ at \$x=100\$", xlabel="t", ylabel="\$\\phi_t\$");
savefig(myplot,"phit_vs_t_at_x=100_Tm3000.png");

temp_phix_vstime_list = phix[startcut:endcut,x_loc];
myplot = plot(t,temp_phix_vstime_list, title = "\$\\phi_x(t)\$ at \$x=100\$", xlabel="t", ylabel="\$\\phi_x\$");
savefig(myplot,"phix_vs_t_at_x=100_Tm3000.png");

myplot = plot(t,g*temp_phit_vstime_list./dt, title = "\$J_x(t)\$ at \$x=100\$", xlabel="t", ylabel="\$\\J_x\$");
savefig(myplot,"current_vs_t_at_x=100_Tm3000.png");

myplot = plot(t,g*tempphi_vstime_list .+ Q, title = "\$E(t)\$ at \$x=100\$", xlabel="t", ylabel="\$E\$");
savefig(myplot,"efield_vs_t_at_x=100_Tm3000.png");


current_at0 = g*phit[:,Int((nv_x-1)/2+1)]./dt;
efield_at0  = g.*phi[:,Int((nv_x-1)/2+1)] .+ Q;
phi_at0     = phi[:,Int((nv_x-1)/2+1)];

myplot = plot(t,phi_at0[startcut:endcut], title = "\$\\phi(t)\$ at \$x=0\$", xlabel="t", ylabel="\$\\phi\$");
savefig(myplot,"phi_vs_t_at_x=0_Tm3000.png");

myplot = plot(t, current_at0[startcut:endcut], title = "\$J_x(t)\$ at \$x=0\$", xlabel="t", ylabel="\$J_x\$");
savefig(myplot,"current_vs_t_at_x=0_Tm3000.png");

myplot = plot(t, efield_at0[startcut:endcut], title = "\$E(t)\$ at \$x=0\$", xlabel="t", ylabel="\$E\$");
savefig(myplot,"efield_vs_t_at_x=0_Tm3000.png");




@save "results_at_x=100_and_x=0.jld2" W g L Tm dx dt Q nv_x nv_t Nv tempphi_vstime_list temp_phit_vstime_list temp_phix_vstime_list current_at0 efield_at0 phi_at0 energycharge energycurrent energyfield
