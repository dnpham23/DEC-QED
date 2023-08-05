using JLD2;
include("mesh.jl"); include("tree.jl"); include("utils.jl"); include("linsolve.jl"); include("amperesolve.jl");include("amperesolve_delta.jl"); include("plotting.jl"); include("computeFields.jl");

xmax =  2000.0;
xmin = -2000.0;
ymax =  2000.0;
ymin = -2000.0;
zmax =  2000.0;
zmin = -2000.0;
Lx_max =  1000.0;
Lx_min = -1000.0;
Ly_max =  1000.0;
Ly_min = -1000.0;
Lz1    =  -400.0;
Lz2    =   0.0;
Lz3    =   400.0;
Nc     = 3;

Sx_max =  400.0;
Sx_min = -400.0;
Sy_max =  400.0;
Sy_min = -400.0;
Sz_max =  400.0;
Sz_min = -400.0;

Nx = 21;           # number of vertices along x
Ny = 21;           # number of vertices along y
Nz = 21;           # number of vertices along z
ne_x = Nx - 1;    # number of x-edge on each row
ne_y = Ny - 1;    # number of y-edge on each column
ne_z = Nz - 1;    # number of z-edge on each column
Ne_x = ne_x*Ny*Nz;   # total number of x-edges
Ne_y = ne_y*Nx*Nz;   # total number of y-edges
Ne_z = ne_z*Nx*Ny;   # total number of z-edges
Ne = Ne_x + Ne_y + Ne_z; # total number of edges
Nv = Nx*Ny*Nz;       # total number of vertices
Nbranch = Nx*ne_y + ne_x + Ne_z;
lx = (xmax-xmin)/ne_x;
ly = (ymax-ymin)/ne_y;
lz = (zmax-zmin)/ne_z;
Nv_xyplane  = Nx*Ny;
Nex_xyplane = ne_x*Ny;
Ney_xyplane = ne_y*Nx;
Ne_xyplane = Nex_xyplane + Ney_xyplane; # doesn't inlcude z edges

# edge info for the current loop
ne_x_loop = Int((Lx_max - Lx_min)/lx);
ne_y_loop = Int((Ly_max - Ly_min)/ly);
ne_per_ring = 2*(ne_x_loop + ne_y_loop);
necurrent = Nc*ne_per_ring;

# edge info for the superconducting piece
ne_scx = Int((Sx_max - Sx_min)/lx);
ne_scy = Int((Sy_max - Sy_min)/ly);
ne_scz = Int((Sz_max - Sz_min)/lz);
Ne_scx = ne_scx*(ne_scy+1)*(ne_scz+1);
Ne_scy = (ne_scx+1)*ne_scy*(ne_scz+1);
Ne_scz = (ne_scx+1)*(ne_scy+1)*ne_scz;
Ne_sc  = Ne_scx + Ne_scy + Ne_scz;
Ne_scx_xyplane = ne_scx*(ne_scy+1);
Ne_scy_xyplane = ne_scy*(ne_scx+1);
Ne_scz_xyplane = (ne_scx+1)*(ne_scy+1);

# edge info for air 
ne_airx_seg = Int((Sx_min - xmin)/lx);
ne_airy_seg = Int((Sy_min - ymin)/ly);
ne_airz_seg = Int((Sz_min - zmin)/lz);

# Some constants
fo       = 1e9;                # 1GHz
to       = 1/fo;               # 1 ns
lo       = 1e-6;               # 1 um
lambda_o = 1e-6;               # 1 um
Qo       = -2*1.60217663e-19;   # 2*electron charge
mo       = 2*9.1093837015e-31; # 2*electron mass
Phi_o    = 2.067833848e-15;    # magnetic flux quantum
eps_o    = 8.8541878128e-12;   # vacuum permitivity
mu_o     = 1.25663706212e-6;   # vacuum magnetic permeability
C        = 299792458; # speed of light

# Input values
Nstep     = 20;
simulatetime = 5e3; # nanoseconds
dt        = simulatetime/(Nstep-1);
alpha     = 0.1; # max current is J = alpha*simulatetime*(Qo/to) = 1.6 uA
eps_air   = 1.0;
eps_sc    = 1.0;
mu_air    = 1.0;
mu_sc     = 1.0;
lambda    = 200.0; # in micrometer
max_iter  = 15;
conv_tol  = 0.05;


# scaling for Ampere's and Gauss's equation
A1       = to^2/(mu_o*eps_o*lo^2);
#A2       = (to*Qo)/(eps_o*lo*Phi_o);
A2       = (to*Qo)/(eps_o*Phi_o); # instead of having lo in denominator, absorb it into Q above
G1       = to^2/(eps_o*mu_o*lambda_o^2);
G2       = Qo*Phi_o*to/(mo*lo^2);
G3       = to*Phi_o/(mu_o*Qo*lo);
S1       = A1;
S2       = to^2/(mu_o*eps_o*lambda_o^2);
S3       = (Qo*to*Phi_o)/(mo*lo^2);
S4       = (to*Qo*Phi_o)/(mo*lo^2);
S5       = (Qo^2*to^2)/(mo*eps_o*lo^3);
J1       = to*Qo/(eps_o*Phi_o*lo); # true J1 scaling
Is1      = (Phi_o*to)/(mu_o*Qo*lo);
Is2      = (eps_o*Phi_o^2)/(mo*lo);
rho1     = mu_o*eps_o*Qo*Phi_o/(mo*to);


# define locations of the current edges
currentloc    = zeros(2,3,necurrent);
# 1st RING
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,m] = [xc Ly_min Lz1; xc+lx Ly_min Lz1]; # currents live on the edges
end
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,ne_x_loop+m] = [xc Ly_max Lz1; xc+lx Ly_max Lz1]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,2*ne_x_loop+m] = [Lx_min yc Lz1; Lx_min yc+ly Lz1]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,2*ne_x_loop+ne_y_loop+m] = [Lx_max yc Lz1; Lx_max yc+ly Lz1]; # currents live on the edges
end
# 2nd RING
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,ne_per_ring+m] = [xc Ly_min Lz2; xc+lx Ly_min Lz2]; # currents live on the edges
end
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,ne_per_ring+ne_x_loop+m] = [xc Ly_max Lz2; xc+lx Ly_max Lz2]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,ne_per_ring+2*ne_x_loop+m] = [Lx_min yc Lz2; Lx_min yc+ly Lz2]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,ne_per_ring+2*ne_x_loop+ne_y_loop+m] = [Lx_max yc Lz2; Lx_max yc+ly Lz2]; # currents live on the edges
end
# 3rd ring
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,2*ne_per_ring+m] = [xc Ly_min Lz3; xc+lx Ly_min Lz3]; # currents live on the edges
end
for m = 1:ne_x_loop
    xc = Lx_min + (m-1)*lx;
    currentloc[:,:,2*ne_per_ring+ne_x_loop+m] = [xc Ly_max Lz3; xc+lx Ly_max Lz3]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,2*ne_per_ring+2*ne_x_loop+m] = [Lx_min yc Lz3; Lx_min yc+ly Lz3]; # currents live on the edges
end
for m = 1:ne_y_loop
    yc = Ly_min + (m-1)*ly;
    currentloc[:,:,2*ne_per_ring+2*ne_x_loop+ne_y_loop+m] = [Lx_max yc Lz3; Lx_max yc+ly Lz3]; # currents live on the edges
end
# define current array
currentamp = ones(necurrent);

e,v  = Mesh3Dcube(xmin,ymin,zmin,Nx,ne_x,Ne_x,Ne_y, Ne_z, Ne,Nv,lx,ly,lz, Nv_xyplane, Nex_xyplane, Ney_xyplane);
vbound, xtol, ytol, ztol  = vboundary3D(v,xmax,xmin,ymax,ymin,zmax,zmin, Nv,Nx,Ny,Nv_xyplane);
v2exmap, v2eymap, v2ezmap = vemap3D(v,Nx,Nv,ne_x,Ne_x, Ne_y, xtol, ytol,ztol,xmin,xmax,ymin,ymax,zmin,zmax, Nex_xyplane, Ney_xyplane, Nv_xyplane);
e_sc = regionsort3D(Nx,ne_x,Ne_x,Ne_y, ne_scx,Ne_scx,Ne_scy,Ne_scz, Ne_sc, Ne_scx_xyplane, Ne_scy_xyplane, Ne_scz_xyplane, Nv_xyplane, Nex_xyplane, Ney_xyplane, ne_airx_seg, ne_airy_seg,ne_airz_seg);
ebound, ebound_all  = eboundary3D(vbound,e,ne_x,ne_y,ne_z,Ne,Nx, Ny, Nz, Ne_x, Ne_y, Nex_xyplane, Ney_xyplane);
invLambda2 = materials(Ne, e_sc, lambda);

tree       = stree3D(ne_x,Nbranch,Ney_xyplane,Ne_x,Ne_y, Ne_z);

# Allocate space to hold the solutions
phixpgrid = zeros(ne_x,Ny,Nz,Nstep);
phiypgrid = zeros(Nx,ne_y,Nz,Nstep);
phizpgrid = zeros(Nx,Ny,ne_z,Nstep);
rhogrid   = zeros(Nx,Ny,Nz,Nstep);
Bxgrid    = zeros(Nx, ne_y, ne_z, Nstep);
Bygrid    = zeros(ne_x, Ny, ne_z, Nstep);
Bzgrid    = zeros(ne_x, ne_y, Nz, Nstep);
Ixgrid    = zeros(ne_x,Ny,Nz,Nstep);
Iygrid    = zeros(Nx,ne_y,Nz,Nstep);
Izgrid    = zeros(Nx,Ny,ne_z,Nstep);

for tn = 1:Nstep
    println("Time step: ", tn, "\n");
    # Field values of the previous timestep
    if tn == 1
        phixp_past1 = zeros(Ne_x); 
        phiyp_past1 = zeros(Ne_y);
        phizp_past1 = zeros(Ne_z);
        phixp_past2 = zeros(Ne_x);
        phiyp_past2 = zeros(Ne_y);
        phizp_past2 = zeros(Ne_z);
        rho_past1   = zeros(Nv);
    elseif tn == 2
        # phip            
        phixp_past1 = reshape(phixpgrid[:,:,:,tn-1],Ne_x);
        phiyp_past1 = reshape(phiypgrid[:,:,:,tn-1],Ne_y);
        phizp_past1 = reshape(phizpgrid[:,:,:,tn-1],Ne_z);
        phixp_past2 = zeros(Ne_x);
        phiyp_past2 = zeros(Ne_y);
        phizp_past2 = zeros(Ne_z);
        # rho
        rho_past1   = reshape(rhogrid[:,:,:,tn-1],Nv);
    else
        # phip
        phixp_past1 = reshape(phixpgrid[:,:,:,tn-1],Ne_x);
        phiyp_past1 = reshape(phiypgrid[:,:,:,tn-1],Ne_y);
        phizp_past1 = reshape(phizpgrid[:,:,:,tn-1],Ne_z);
        phixp_past2 = reshape(phixpgrid[:,:,:,tn-2],Ne_x);
        phiyp_past2 = reshape(phiypgrid[:,:,:,tn-2],Ne_y);
        phizp_past2 = reshape(phizpgrid[:,:,:,tn-2],Ne_z);
        # rho
        rho_past1   = reshape(rhogrid[:,:,:,tn-1],Nv);
    end
    # assign current values
    I = alpha*(tn-1)*dt;
    # 1st ring
    currentamp[1:ne_x_loop] = I*ones(ne_x_loop);
    currentamp[ne_x_loop+1:2*ne_x_loop] = -1.0*I*ones(ne_x_loop);
    currentamp[2*ne_x_loop+1:2*ne_x_loop+ne_y_loop] = -1.0*I*ones(ne_y_loop);
    currentamp[2*ne_x_loop+ne_y_loop+1:ne_per_ring] = I*ones(ne_y_loop);
    # 2nd ring
    currentamp[ne_per_ring+1 : ne_per_ring+ne_x_loop] = I*ones(ne_x_loop);
    currentamp[ne_per_ring+ne_x_loop+1 : ne_per_ring+2*ne_x_loop] = -1.0*I*ones(ne_x_loop);
    currentamp[ne_per_ring+2*ne_x_loop+1 : ne_per_ring+2*ne_x_loop+ne_y_loop] = -1.0*I*ones(ne_y_loop);
    currentamp[ne_per_ring+2*ne_x_loop+ne_y_loop+1 : 2*ne_per_ring] = I*ones(ne_y_loop);
    # 3rd ring
    currentamp[2*ne_per_ring+1 : 2*ne_per_ring+ne_x_loop] = I*ones(ne_x_loop);
    currentamp[2*ne_per_ring+ne_x_loop+1 : 2*ne_per_ring+2*ne_x_loop] = -1.0*I*ones(ne_x_loop);
    currentamp[2*ne_per_ring+2*ne_x_loop+1 : 2*ne_per_ring+2*ne_x_loop+ne_y_loop] = -1.0*I*ones(ne_y_loop);
    currentamp[2*ne_per_ring+2*ne_x_loop+ne_y_loop+1 : 3*ne_per_ring] = I*ones(ne_y_loop);
    # assign current values to the correct edges
    ecurrent   = current3D(v,currentloc,currentamp,xtol,ytol,ztol, Nv, Nx, Ne,e);

    # solve for rho
    rho_current = chargesolve_timedomain(G2, G3, dt,lx,ly,lz,phixp_past1, phiyp_past1, phizp_past1, v2exmap,v2eymap,v2ezmap, Nv, rho_past1, invLambda2,e_sc,e, Ne_x, Ne_y);
    rhogrid[:,:,:,tn] = reshape(rho_current,Nx,Ny,Nz);
    # construct the matrices
    AMat = Ampere_AMat3D(S4,dt,Ne,ne_x, Nx,Ne_x,Ne_y,ebound, e_sc, lx,ly,lz,Ny,ne_y, Nex_xyplane, Ney_xyplane, Nv_xyplane);
    BMat = Ampere_BMat(dt,Ne,Nx, Ny,ne_x,ebound_all);
    Ccol = Ampere_Ccol3D(J1, S1, S2, S4, S5, dt,ne_x, ne_y, Nx, Ne_x, Ne_y, ebound_all, phixp_past1, phiyp_past1, phizp_past1, phixp_past2, phiyp_past2, phizp_past2, rho_past1, e_sc, Ny, lx, ly, lz, invLambda2,Ne, ecurrent, Nex_xyplane, Ney_xyplane, Nv_xyplane, e); 

    # solve for phi_p
    MMat      = Ampere_MMat(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1, ebound_all, Ne);
    Kcol      = Ampere_Kcol(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1);
    phip_next = ampere_timedomain_delta3D(MMat,Kcol,phixp_past1, phiyp_past1, phizp_past1);

    # Need to sort phip_next into phip_x and phip_y
    phip_next_x = phip_next[1:Ne_x];
    phip_next_y = phip_next[Ne_x+1:Ne_x+Ne_y];
    phip_next_z = phip_next[Ne_x+Ne_y+1:end];    

    Bx, By, Bz      = computeBfield3D(phix,phiy,phiz,Nx, Ny, Nz,ne_x, ne_y, ne_z, lx, ly, lz,Nex_xyplane, Ney_xyplane, Nv_xyplane);
    phixpgrid[:,:,:,tn] = reshape(phip_next_x,ne_x,Ny,Nz);
    phiypgrid[:,:,:,tn] = reshape(phip_next_y,Nx,ne_y,Nz);
    phizpgrid[:,:,:,tn] = reshape(phip_next_z,Nx,Ny,ne_z);
    Bxgrid[:,:,:,tn]    = reshape(Bx, Nx, ne_y, ne_z);
    Bygrid[:,:,:,tn]    = reshape(By, ne_x, Ny, ne_z);
    Bzgrid[:,:,:,tn]    = reshape(Bz, ne_x, ne_y, Nz);
    rhogrid[:,:,:,tn]   = reshape(rho,Nx,Ny,Nz);
end

@save "SimulationData3.jld2" phixpgrid phiypgrid phizpgrid Bx By Bz rhogrid
