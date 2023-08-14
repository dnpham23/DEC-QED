using JLD2;
#using IterativeSolvers;
include("mesh.jl"); include("tree.jl"); include("utils.jl"); include("linsolve.jl"); include("amperesolve.jl");include("amperesolve_delta.jl"); include("computeFields.jl");

# dimensions of the problem (all lengths are scaled by penetration depth lambda_L)
xmax =  24.0;
xmin = -24.0;
ymax =  24.0;
ymin = -24.0;
zmax =  24.0;
zmin = -24.0;
xmax_cav =  18.0; 
xmin_cav = -18.0; 
ymax_cav =  18.0;
ymin_cav = -18.0; 
zmax_cav =  18.0;
zmin_cav = -18.0;

Cloc1 = -4.0;
Cloc2 =  4.0;

Nx = 25;           # number of vertices along x
Ny = 25;           # number of vertices along y
Nz = 25;           # number of vertices along z

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

# info about sc edges
ne_airx = Int((xmax_cav-xmin_cav)/lx);
ne_airy = Int((ymax_cav-ymin_cav)/ly);
ne_airz = Int((zmax_cav-zmin_cav)/lz);
Ne_airx = ne_airx*(ne_airy+1)*(ne_airz+1); # number of x-edges in air
Ne_airy = ne_airy*(ne_airx+1)*(ne_airz+1); # number of y-edges in air
Ne_airz = (ne_airx+1)*(ne_airy+1)*ne_airz; # number of z-edges in air
Ne_air  = Ne_airx + Ne_airy + Ne_airz;
Nv_air  = (ne_airx+1)*(ne_airy+1)*(ne_airz+1);

Ne_airx_xyplane = ne_airx*(ne_airy+1);
Ne_airy_xyplane = ne_airy*(ne_airx+1);
Ne_airz_xyplane = (ne_airx+1)*(ne_airy+1);
Ne_airx_xzplane = ne_airx*(ne_airz+1);
Ne_airz_xzplane = ne_airz*(ne_airx+1);
Ne_airy_yzplane = ne_airy*(ne_airz+1);
Ne_airz_yzplane = ne_airz*(ne_airy+1);

ne_scx_seg  = Int((xmax-xmax_cav)/lx);
ne_scy_seg  = Int((ymax-ymax_cav)/ly);
ne_scz_seg  = Int((zmax-zmax_cav)/lz);
Ne_scx      = 2*ne_x*Ny*(ne_scz_seg+1) + 2*ne_scx_seg*Ny*(ne_airz-1) + 2*ne_airx*(ne_scy_seg+1)*(ne_airz-1); # number of x-edges in the superconductor
Ne_scy      = 2*ne_y*Nx*(ne_scz_seg+1) + 2*ne_scy_seg*Nx*(ne_airz-1) + 2*ne_airy*(ne_scx_seg+1)*(ne_airz-1); # number of y-edges in the superconductor
Ne_scz      = 2*Nx*Ny*ne_scz_seg + 2*Nx*(ne_scy_seg+1)*ne_airz + 2*(ne_scx_seg+1)*(ne_airy-1)*ne_airz; # number of z-edges in the superconductor
Ne_sc  = Ne_scx + Ne_scy + Ne_scz;
#Nv_sc  = Nv_sc_xyplane*(ne_scz+1);


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

# Input values
Q         = 1.0*1e3; 
#Q         = 1.0e8;
freq      = 1e4;
Period    = 1/freq; # 
Omega     = 2*pi*freq;
Nstep     = 201;
nperiod   = 5;
simulatetime = nperiod*Period;
dt        = simulatetime/(Nstep-1);
lambda    = 1.0; # in micrometer
max_iter  = 15;
reltol    = 1e-4;

# define locations of charges and currents
chargeloc     = [0 0 Cloc1;0 0 Cloc2];
chargeamp     = [0.0; 0.0];
currentendpts = chargeloc;
Ncurrent      = Int(div(currentendpts[2,3]-currentendpts[1,3],ly));
currentloc    = zeros(2,3,Ncurrent);
for m = 1:Ncurrent
    zc1 = currentendpts[1,3] + (m-1)*lz;
    zc2 = currentendpts[1,3] + m*lz;
    currentloc[:,:,m] = [0.0 0.0 zc1;0.0 0.0 zc2]; # currents live on the edges
end
# define current array
#currentamp = ones(Ncurrent);

# define the mesh
e,v  = Mesh3Dcube(xmin,ymin,zmin,Nx,ne_x,Ne_x,Ne_y, Ne_z, Ne,Nv,lx,ly,lz, Nv_xyplane, Nex_xyplane, Ney_xyplane);
# find boundary nodes
vbound, xtol, ytol, ztol  = vboundary3D(v,xmax,xmin,ymax,ymin,zmax,zmin, Nv,Nx,Ny,Nv_xyplane);
# mapping v->e and e->v
v2exmap, v2eymap, v2ezmap = vemap3D(v,Nx,Nv,ne_x,Ne_x, Ne_y, xtol, ytol,ztol,xmin,xmax,ymin,ymax,zmin,zmax, Nex_xyplane, Ney_xyplane, Nv_xyplane);
# find the edges in the sc region
e_sc, e_air = regionsort3D_cav(e,v,Nx,Ny,Nv,ne_x,ne_y,Ne_x,Ne_y, lx,ly, ne_airx, ne_airy, ne_airz, Ne_airx, Ne_airy,  Ne_airz, Ne_air, Ne_airx_xyplane, Ne_airy_xyplane, Ne_airz_xyplane, Ne_airx_xzplane, Ne_airz_xzplane, Ne_airy_yzplane, Ne_airz_yzplane, Ne_scx, Ne_scy, Ne_scz, Ne_sc, Nv_xyplane, Nex_xyplane, Ney_xyplane, ne_scx_seg, ne_scy_seg,ne_scz_seg);
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

BLAS.set_num_threads(1);
for tn = 1:Nstep
#for tn = 1:5
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

    currentamp = -Q*Omega*cos(2*Omega*(tn-1)*dt)*ones(Ncurrent);    
    # assign current values to the correct edges
    ecurrent   = current3D(v,currentloc,currentamp,xtol,ytol,ztol, Nv, Nx, Ne,e);

    # solve for rho
    rho_current = chargesolve_timedomain_cav(G2, G3, dt,lx,ly,lz,phixp_past1, phiyp_past1, phizp_past1, v2exmap,v2eymap,v2ezmap, Nv, rho_past1, invLambda2,e_sc,e, Ne_x, Ne_y, ecurrent);
    rhogrid[:,:,:,tn] = reshape(rho_current,Nx,Ny,Nz);
    # construct the matrices
    AMat = Ampere_AMat3D(S4,dt,Ne,ne_x, Nx,Ne_x,Ne_y,ebound, e_sc, lx,ly,lz,Ny,ne_y, Nex_xyplane, Ney_xyplane, Nv_xyplane);
    BMat = Ampere_BMat(dt,Ne,Nx, Ny,ne_x,ebound_all);
    Ccol = Ampere_Ccol3D(J1, S1, S2, S4, S5, dt,ne_x, ne_y, Nx, Ne_x, Ne_y, ebound_all, phixp_past1, phiyp_past1, phizp_past1, phixp_past2, phiyp_past2, phizp_past2, rho_past1, e_sc, Ny, lx, ly, lz, invLambda2,Ne, ecurrent, Nex_xyplane, Ney_xyplane, Nv_xyplane, e); 

    # solve for phi_p
    MMat      = Ampere_MMat(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1, ebound_all, Ne);
    Kcol      = Ampere_Kcol(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1);
    ###phip_next = ampere_timedomain_delta3D(MMat,Kcol,phixp_past1, phiyp_past1, phizp_past1);
    phip_next = ampere_timedomain_delta3D_gmres(MMat,Kcol,phixp_past1, phiyp_past1, phizp_past1, reltol) 

    # Need to sort phip_next into phip_x and phip_y
    phip_next_x = phip_next[1:Ne_x];
    phip_next_y = phip_next[Ne_x+1:Ne_x+Ne_y];
    phip_next_z = phip_next[Ne_x+Ne_y+1:end];    

    Bx, By, Bz      = computeBfield3D(phip_next_x, phip_next_y, phip_next_z, Nx, Ny, Nz,ne_x, ne_y, ne_z, lx, ly, lz, Nex_xyplane, Ney_xyplane, Nv_xyplane);

    phixpgrid[:,:,:,tn] = reshape(phip_next_x,ne_x,Ny,Nz);
    phiypgrid[:,:,:,tn] = reshape(phip_next_y,Nx,ne_y,Nz);
    phizpgrid[:,:,:,tn] = reshape(phip_next_z,Nx,Ny,ne_z);
    Bxgrid[:,:,:,tn]    = reshape(Bx, Nx, ne_y, ne_z);
    Bygrid[:,:,:,tn]    = reshape(By, ne_x, Ny, ne_z);
    Bzgrid[:,:,:,tn]    = reshape(Bz, ne_x, ne_y, Nz);
    rhogrid[:,:,:,tn]   = reshape(rho_current,Nx,Ny,Nz);
end

@save "SimulationData3.jld2" phixpgrid phiypgrid phizpgrid Bxgrid Bygrid Bzgrid rhogrid Nstep lx ly lz
