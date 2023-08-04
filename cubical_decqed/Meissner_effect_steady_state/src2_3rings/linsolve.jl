function reducedlinsolve(Y::Array{Float64,2},x::Array{Float64,2},vcharge::Array{Float64,2},ebound_all::Array{Int,1},tree::Array{Int,1})
    Mat = Y*x;
    for i = 1:length(tree)
        if (tree[i] in ebound_all)
            Mat[i,:]  .= 0.0;
            Mat[i,i]   = 1.0;
            vcharge[i] = 0.0;
        end
    end
    psireduced = inv(Mat)*(-vcharge);
    return psireduced;
end

function interpPsi(psireduced::Array{Float64,2},x::Array{Float64,2},Ne_x::Int64)
    psi  = x*psireduced;
    psix = vec(psi[1:Ne_x,:]);
    psiy = vec(psi[Ne_x+1:end,:]);
    return psi, psix, psiy;
end

function linsolve_ssPhi(HMat::Array{Float64,2},z::Array{Float64,2},Kcol::Array{Float64,1}, ebound::Array{Int,1}, Ne_x::Int64,Is1::Float64, invLambda2::Array{Float64,1}, J::Float64, ne_x::Int64)
    
    RMat         = z'*HMat*z;
    Kcol_reduced = z'*Kcol; 
    
    ## set x-edge at the corners to 0
    #rc1 = (ne_y-1)*ne_x + 1;
    #rc2 = ne_y*ne_x;
    #RMat[rc1,:]   .= 0.0;
    #RMat[rc1,rc1]  = 1.0;
    #RMat[rc2,:]   .= 0.0;
    #RMat[rc2,rc2]  = 1.0;
    #Kcol_reduced[rc1] = 0.0;
    #Kcol_reduced[rc2] = 0.0;
    
    ## set boundary condition
    ## first the x-edges at the vertical boundaries
    #for i = 2:(ne_y+1)
    #	xind1 = (i-1)*ne_x + 1;
    #	xind2 = i*ne_x;
    #	if (!(xind1 in pedge_x))&&(!(xind2 in pedge_x))
    #	    rind1 = xind1 - ne_x;
    #	    rind2 = xind2 - ne_x;
    #	    RMat[rind1,:] .= 0.0;
    #	    RMat[rind2,:] .= 0.0;
    #	    RMat[rind1, rind1] = 1.0;
    #	    RMat[rind2, rind2] = 1.0;
    #	    Kcol_reduced[rind1] = 0.0;
    #	    Kcol_reduced[rind2] = 0.0;
    #	end
    #end
    ## then the x-edges at the upper horizontal boundary
    #for i = 1:ne_x
    #	rind = (ne_y-1)*ne_x + i;
    #	RMat[rind,:]    .= 0.0;
    #	RMat[rind, rind] = 1.0;
    #	Kcol_reduced[rind] = 0.0;
    #end

    phireduced   = RMat\Kcol_reduced;
    phip_ss      = z*phireduced;
    phip_ss_x    = phip_ss[1:Ne_x];
    phip_ss_y    = phip_ss[Ne_x+1:Ne_x+Ne_y];
    phip_ss_z    = phip_ss[Ne_x+Ne_y+1:end];
    
    return phip_ss_x, phip_ss_y, phip_ss_z;
end

function linsolve_ssPhi3D(RMat::Array{Float64,2},z::Array{Float64,2},Kcol_reduced::Array{Float64,1}, ebound::Array{Int,1}, Ne_x::Int64,Is1::Float64, invLambda2::Array{Float64,1}, J::Float64, ne_x::Int64)

    phireduced   = RMat\Kcol_reduced;
    phip_ss      = z*phireduced;
    phip_ss_x    = phip_ss[1:Ne_x];
    phip_ss_y    = phip_ss[Ne_x+1:Ne_x+Ne_y];
    phip_ss_z    = phip_ss[Ne_x+Ne_y+1:end];
    
    return phip_ss_x, phip_ss_y, phip_ss_z;
end

function linsolve_ssPhi3D_gmres(RMat::Array{Float64,2},z::Array{Float64,2},Kcol_reduced::Array{Float64,1}, ebound::Array{Int,1}, Ne_x::Int64,Is1::Float64, invLambda2::Array{Float64,1}, J::Float64, ne_x::Int64, reltol::Float64)
    
    phireduced, history = gmres(RMat, Kcol_reduced, log=true, restart=10, maxiter=20, reltol=reltol);
    phip_ss      = z*phireduced;
    phip_ss_x    = phip_ss[1:Ne_x];
    phip_ss_y    = phip_ss[Ne_x+1:Ne_x+Ne_y];
    phip_ss_z    = phip_ss[Ne_x+Ne_y+1:end];
    
    return phip_ss_x, phip_ss_y, phip_ss_z, history;
end
