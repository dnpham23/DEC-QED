#function main()
    using JLD2;
    #using IterativeSolvers;
    include("mesh.jl"); include("tree.jl"); include("utils.jl"); include("linsolve.jl"); include("amperesolve.jl");include("amperesolve_delta.jl"); include("plotting.jl"); include("computeFields.jl");

# dimensions of the problem (all lengths are scaled by penetration depth lambda_L)
xmax =  30.0;
xmin = -30.0;
ymax =  30.0;
ymin = -30.0;
zmax =  16.0;
zmin = -16.0;
xleft_scring_inner = -14.0;
xleft_scring_outer = -24.0;
xright_scring_inner = 14.0;    
xright_scring_outer = 24.0;  
ybottom_scring_inner = -14.0;
ybottom_scring_outer = -24.0;
ytop_scring_inner    =  14.0;
ytop_scring_outer    =  24.0;
zmax_scring          =  6.0;
zmin_scring          = -6.0;

Lx_max =  8.0;
Lx_min = -8.0;
Ly_max =  8.0;
Ly_min = -8.0;
Lz_max =  6.0;
Lz_min = -6.0;

Nx = 31;           # number of vertices along x
Ny = 31;           # number of vertices along y
Nz = 17;           # number of vertices along z


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
ne_scx_seg  = Int((xleft_scring_inner-xleft_scring_outer)/lx);
ne_scy_seg  = Int((ybottom_scring_inner-ybottom_scring_outer)/ly);
ne_scz      = Int((zmax_scring-zmin_scring)/lz);
ne_scx_len  = Int((xright_scring_outer-xleft_scring_outer)/lx);
ne_scy_len  = Int((ytop_scring_outer-ybottom_scring_outer)/ly);
ne_scx_inlen = Int((xright_scring_inner-xleft_scring_inner)/lx);
ne_scy_inlen = Int((ytop_scring_inner-ybottom_scring_inner)/ly);
ne_scx_horblock = ne_scx_len*(ne_scy_seg+1);   # horblock means horizontal in xy plane
ne_scx_verblock = ne_scx_seg*(ne_scy_inlen-1); # verblock means vertical in xy plane
ne_scy_verblock = ne_scy_len*(ne_scx_seg+1);   # verblock means vertical in xy plane
ne_scy_horblock = ne_scy_seg*(ne_scx_inlen-1); # horblock means horizontal in xy plane
nv_sc_horblock  = (ne_scx_len+1)*(ne_scy_seg+1);
nv_sc_verblock  = (ne_scx_seg+1)*(ne_scy_inlen-1);
nv_scbound      = 2*(ne_scx_len+ne_scy_len+ne_scx_inlen+ne_scy_inlen)*(ne_scz+1);

Ne_jj  = (ne_scx_seg+1)*(ne_scz+1);
Ne_scx_xyplane = 2*(ne_scy_seg+1)*ne_scx_len + 2*ne_scx_seg*(ne_scy_inlen-1); # number of x-edges in sc on one xy plane
Ne_scy_xyplane = 2*(ne_scx_seg+1)*ne_scy_len + 2*ne_scy_seg*(ne_scx_inlen-1); # number of y-edges in sc on one xy plane
Nv_sc_xyplane  = 2*(ne_scy_seg+1)*(ne_scx_len+1) + 2*(ne_scx_seg+1)*(ne_scy_inlen-1); 
Ne_scx = Ne_scx_xyplane*(ne_scz+1); # number of x-edges in the superconductor
Ne_scy = Ne_scy_xyplane*(ne_scz+1); # number of y-edges in the superconductor
Ne_scz = Nv_sc_xyplane*ne_scz; # number of z-edges in the superconductor
Ne_sc  = Ne_scx + Ne_scy + Ne_scz;
Nv_sc  = Nv_sc_xyplane*(ne_scz+1);

ne_airx_seg = Int((xleft_scring_outer - xmin)/lx);
ne_airy_seg = Int((ybottom_scring_outer - ymin)/ly);;
ne_airz_seg = Int((zmin_scring - zmin)/lz);

# edge info for the current loop
ne_x_loop  = Int((Lx_max - Lx_min)/lx);
ne_y_loop  = Int((Ly_max - Ly_min)/ly);
nz_loop    = Int((Lz_max - Lz_min)/lz)+1;
ne_perloop = 2*(ne_x_loop + ne_y_loop);
necurrent  = ne_perloop*nz_loop;


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
Nstep     = 101;
simulatetime = 3.5e-4;
dt        = simulatetime/(Nstep-1); 
alpha     = 2.865677950146584e8;
Jc        = 10.0*alpha*simulatetime; # critical current per z length
lambda    = 10.0;
conv_tol  = 0.05;
reltol    = 1e-4;

# define locations of the current edges
currentloc    = zeros(2,3,necurrent);
for i = 1:nz_loop
    zc = Lz_min + (i-1)*lz;
    for m = 1:ne_x_loop
        xc = Lx_min + (m-1)*lx;
        currentloc[:,:,(i-1)*ne_perloop + m] = [xc Ly_min zc; xc+lx Ly_min zc]; # currents live on the edges
    end
    for m = 1:ne_x_loop
        xc = Lx_min + (m-1)*lx;
        currentloc[:,:,(i-1)*ne_perloop + ne_x_loop+m] = [xc Ly_max zc; xc+lx Ly_max zc]; # currents live on the edges
    end
    for m = 1:ne_y_loop
        yc = Ly_min + (m-1)*ly;
        currentloc[:,:,(i-1)*ne_perloop + 2*ne_x_loop+m] = [Lx_min yc zc; Lx_min yc+ly zc]; # currents live on the edges
    end
    for m = 1:ne_y_loop
        yc = Ly_min + (m-1)*ly;
        currentloc[:,:,(i-1)*ne_perloop + 2*ne_x_loop+ne_y_loop+m] = [Lx_max yc zc; Lx_max yc+ly zc]; # currents live on the edges
    end
end
# define current array
currentamp = ones(necurrent);

# define the JJ edges
e_jj = Array{Int,1}(undef, Ne_jj);
for i = 1:Ne_jj
    rownum_jj = div(i-1, ne_scx_seg+1)+1;
    colnum_jj = (i-1)%(ne_scx_seg+1) + 1;
    e_jj[i] = Ne_x + (rownum_jj - 1 + ne_airz_seg)*Ney_xyplane + (ne_airy_seg + Int(floor(ne_scy_len/2)) -1)*Nx + ne_airx_seg + colnum_jj;
end

# define the mesh
e,v  = Mesh3Dcube(xmax,xmin,ymax,ymin,zmax,zmin,Nx,Ny,Nz, ne_x,ne_y,ne_z,Ne_x,Ne_y, Ne_z, Ne,Nv,lx,ly,lz, Nv_xyplane, Nex_xyplane, Ney_xyplane);
# find boundary nodes
vbound, xtol, ytol, ztol  = vboundary3D(v,xmax,xmin,ymax,ymin,zmax,zmin, Nv,Nx,Ny,Nv_xyplane);
# mapping v->e and e->v
v2exmap, v2eymap, v2ezmap = vemap3D(e,v,vbound,Nx,Ny,Nz,Nv,ne_x,Ne_x, Ne_y, xtol, ytol,ztol,xmin,xmax, ymin,ymax, zmin,zmax, Nex_xyplane, Ney_xyplane, Nv_xyplane);
# find the edges in the sc region
e_sc = regionsort3D_scring(e,v,Nx,Ny,Nv,ne_x,ne_y,Ne_x,Ne_y, lx,ly,Ne_scx, Ne_scy,Ne_scz, Ne_sc, Nv_xyplane, Nex_xyplane, Ney_xyplane, ne_airx_seg, ne_airy_seg,ne_airz_seg, ne_scx_horblock, ne_scx_verblock, ne_scy_verblock, ne_scy_horblock, nv_sc_horblock, nv_sc_verblock, Ne_scx_xyplane, Ne_scy_xyplane, Nv_sc_xyplane);
ebound, ebound_all  = eboundary3D(vbound,e,ne_x,ne_y,ne_z,Ne,Nx, Ny, Nz, Ne_x, Ne_y, Nex_xyplane, Ney_xyplane);
invLambda2 = materials_scringJJ(Ne, e_sc, e_jj, lambda);

tree       = stree3D(e,Nx,ne_x,ne_y,Nbranch,Ney_xyplane,Ne_x,Ne_y,Ne_z);
#x          = reducedMat(tree,Ne,Nbranch,ne_x,ne_y,Ne_x,Nx); # NEEDS 3D VERSION
#Y          = GaussMat(v2exmap,v2eymap,Nv,Ne,Nx,lx,ly);      # NEEDS 3D VERSION

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

    J = alpha*(tn-1)*dt; # current density (current per cross-section area)
    for i = 1:nz_loop
        base_ne = (i-1)*ne_perloop; 
        currentamp[base_ne+1:base_ne+ne_x_loop] = J*ones(ne_x_loop);
        currentamp[base_ne+ne_x_loop+1:base_ne+2*ne_x_loop] = -1.0*J*ones(ne_x_loop);
        currentamp[base_ne+2*ne_x_loop+1:base_ne+2*ne_x_loop+ne_y_loop] = -1.0*J*ones(ne_y_loop);
        currentamp[base_ne+2*ne_x_loop+ne_y_loop+1:base_ne+2*ne_x_loop+2*ne_y_loop] = J*ones(ne_y_loop);
    end
    
    
    # assign current values to the correct edges
    ecurrent   = current3D(v,currentloc,currentamp,xtol,ytol,ztol, Nv, Nx, Ne,e);
    ecurrentJJ = currentJJ(e_jj,phiyp_past1,Ne_x, Ne_y, e, Jc);

    # solve for rho
    rho_current = chargesolve_timedomain(G2, G3, dt,lx,ly,lz,phixp_past1, phiyp_past1, phizp_past1, v2exmap,v2eymap,v2ezmap, Nv, rho_past1, invLambda2,e_sc,e, Ne_x, Ne_y);
    rhogrid[:,:,:,tn] = reshape(rho_current,Nx,Ny,Nz);
    # construct the matrices
    AMat = Ampere_AMat3D(S4,dt,Ne,ne_x, Nx,Ne_x,Ne_y,ebound, e_sc, lx,ly,lz,Ny,ne_y, Nex_xyplane, Ney_xyplane, Nv_xyplane);
    BMat = Ampere_BMat(dt,Ne,Nx, Ny,ne_x,ebound_all);
    Ccol = Ampere_Ccol3D_2(J1, S1, S2, S4, S5, dt,ne_x, ne_y, Nx, Ne_x, Ne_y, ebound_all, phixp_past1, phiyp_past1, phizp_past1, phixp_past2, phiyp_past2, phizp_past2, rho_past1, e_sc, Ny, lx, ly, lz, invLambda2,Ne, ecurrent, Nex_xyplane, Ney_xyplane, Nv_xyplane, e, ecurrentJJ); 

    # solve for phi_p
    MMat      = Ampere_MMat(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1, ebound_all, Ne);
    Kcol      = Ampere_Kcol(AMat,BMat,Ccol,phixp_past1, phiyp_past1, phizp_past1);
    ###phip_next = ampere_timedomain_delta3D(MMat,Kcol,phixp_past1, phiyp_past1, phizp_past1);
    phip_next = ampere_timedomain_delta3D_gmres(MMat,Kcol,phixp_past1, phiyp_past1, phizp_past1, reltol) 

    # Need to sort phip_next into phip_x and phip_y
    phip_next_x = phip_next[1:Ne_x];
    phip_next_y = phip_next[Ne_x+1:Ne_x+Ne_y];
    phip_next_z = phip_next[Ne_x+Ne_y+1:end];    

    ##phix, phiy  = computePhi(phix_past1, phiy_past1, phip_next_x, phip_next_y, phixp_past1, phiyp_past1, psixpast1, psiypast1, S4,dt,ne_x, Nx, Ne_x,Ne_y, ebound_all,lx,ly);
    Bx, By, Bz      = computeBfield3D(phip_next_x, phip_next_y, phip_next_z, Ne, Nx, Ny, Nz,ne_x, ne_y, ne_z, lx, ly, lz,Nex_xyplane, Ney_xyplane, Nv_xyplane);
    ##Ix, Iy = computeCurrent(phip_next_x,phip_next_y,psix,psiy, ecurrent, lx, ly, eps_list, mu_list, invLambda2, Is1, Is2,ebound_all,Ne_x, Ne_y, Nx, ne_x,lo);

    phixpgrid[:,:,:,tn] = reshape(phip_next_x,ne_x,Ny,Nz);
    phiypgrid[:,:,:,tn] = reshape(phip_next_y,Nx,ne_y,Nz);
    phizpgrid[:,:,:,tn] = reshape(phip_next_z,Nx,Ny,ne_z);
    Bxgrid[:,:,:,tn]    = reshape(Bx, Nx, ne_y, ne_z);
    Bygrid[:,:,:,tn]    = reshape(By, ne_x, Ny, ne_z);
    Bzgrid[:,:,:,tn]    = reshape(Bz, ne_x, ne_y, Nz);
    #Ixgrid[:,:,tn]    = reshape(Ix,ne_x,Ny);
    #Iygrid[:,:,tn]    = reshape(Iy,Nx,ne_y);
    rhogrid[:,:,:,tn]   = reshape(rho_current,Nx,Ny,Nz);
end

#return Lk, Ax_ss_grid, Ay_ss_grid, Ax_ss_vec, Ay_ss_vec;
@save "SimulationData3.jld2" phixpgrid phiypgrid phizpgrid Bxgrid Bygrid Bzgrid rhogrid Nstep lx ly lz
#end
