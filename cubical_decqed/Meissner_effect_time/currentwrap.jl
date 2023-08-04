for tn = 1:Nstep*nperiod
    psix = reshape(psixgrid[:,:,tn],Ne_x);
    psiy = reshape(psiygrid[:,:,tn-1],Ne_y); 
    phixp = reshape(phixpgrid[:,:,tn],Ne_x);
    phiyp = reshape(phiypgrid[:,:,tn],Ne_y);
    Ixgrid = zeros(ne_x,Ny,Nstep*nperiod);
    Iygrid = zeros(Nx,ne_y,Nstep*nperiod);
    Ix, Iy = computeCurrent(phixp,phiyp,psix,psiy, ecurrent, lx, ly, eps_list, mu_list, invLambda2, Is1, Is2,ebound_all,Ne_x, Ne_y, Nx, ne_x);
    Ixgrid[:,:,tn] = reshape(Ix,ne_x,Ny);
    Iygrid[:,:,tn] = reshape(Iy,Nx,ne_y);
end