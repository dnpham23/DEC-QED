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

function chargesolve_timedomain(G2::Float64, G3::Float64, dt::Float64,lx::Float64,ly::Float64,lz::Float64,phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, phizp_past1::Array{Float64,1}, v2exmap::Array{Int,2},v2eymap::Array{Int,2},v2ezmap::Array{Int,2}, Nv::Int64, rho_past1::Array{Float64,1},invLambda2::Array{Float64,1},e_sc::Array{Int,1},e::Array{Int,2}, Ne_x::Int64, Ne_y::Int64)
    rho = zeros(Nv);
    for i = 1:Nv
        xminus_edge = v2exmap[i,1];
        xplus_edge  = v2exmap[i,2];
        yminus_edge = v2eymap[i,1];
        yplus_edge  = v2eymap[i,2];
        zminus_edge = v2ezmap[i,1];
        zplus_edge  = v2ezmap[i,2];
        
        temp = dt/(lx*ly*lz);
        rho_temp = 0.0;
        if (xminus_edge!=0)
            if (xminus_edge in e_sc)
                v1 = e[xminus_edge,1];
                v2 = e[xminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += -(G3*invLambda2[xminus_edge] + G2*rho_avg)*(ly*lz/lx)*phixp_past1[xminus_edge];
            end
        end
        
        if (xplus_edge!=0)
            if (xplus_edge in e_sc)
                v1 = e[xplus_edge,1];
                v2 = e[xplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += (G3*invLambda2[xplus_edge] + G2*rho_avg)*(ly*lz/lx)*phixp_past1[xplus_edge];
            end
        end
        
        if (yminus_edge != 0)
            if (yminus_edge in e_sc)
                v1 = e[yminus_edge,1];
                v2 = e[yminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += -(G3*invLambda2[yminus_edge] + G2*rho_avg)*(lx*lz/ly)*phiyp_past1[yminus_edge-Ne_x];
            end
        end
        
        if (yplus_edge != 0)
            if (yplus_edge in e_sc)
                v1 = e[yplus_edge,1];
                v2 = e[yplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += (G3*invLambda2[yplus_edge] + G2*rho_avg)*(lx*lz/ly)*phiyp_past1[yplus_edge-Ne_x];
            end
        end

        if (zminus_edge != 0)
            if (zminus_edge in e_sc)
                v1 = e[zminus_edge,1];
                v2 = e[zminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += -(G3*invLambda2[zminus_edge] + G2*rho_avg)*(lx*ly/lz)*phizp_past1[zminus_edge-Ne_x-Ne_y];
            end
        end
        
        if (zplus_edge != 0)
            if (zplus_edge in e_sc)
                v1 = e[zplus_edge,1];
                v2 = e[zplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp += (G3*invLambda2[zplus_edge] + G2*rho_avg)*(lx*ly/lz)*phizp_past1[zplus_edge-Ne_x-Ne_y];
            end
        end
        
        rho[i] = rho_temp*temp + rho_past1[i];
    end
    
    return rho;    
end

function chargesolve_timedomain_cav(G2::Float64, G3::Float64, dt::Float64,lx::Float64,ly::Float64,lz::Float64,phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, phizp_past1::Array{Float64,1}, v2exmap::Array{Int,2},v2eymap::Array{Int,2},v2ezmap::Array{Int,2}, Nv::Int64, rho_past1::Array{Float64,1},invLambda2::Array{Float64,1},e_sc::Array{Int,1},e::Array{Int,2}, Ne_x::Int64, Ne_y::Int64, ecurrent::Array{Float64,2})
    rho = zeros(Nv);
    for i = 1:Nv
        xminus_edge = v2exmap[i,1];
        xplus_edge  = v2exmap[i,2];
        yminus_edge = v2eymap[i,1];
        yplus_edge  = v2eymap[i,2];
        zminus_edge = v2ezmap[i,1];
        zplus_edge  = v2ezmap[i,2];
        
        temp = dt/(lx*ly*lz);
        rho_temp  = 0.0;
        divJ_temp = 0.0; 
        if (xminus_edge!=0)
            if (xminus_edge in e_sc)
                v1 = e[xminus_edge,1];
                v2 = e[xminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += -(G3*invLambda2[xminus_edge] + G2*rho_avg)*(ly*lz/lx)*phixp_past1[xminus_edge];
                divJ_temp += ecurrent[xminus_edge];
            end
        end
        
        if (xplus_edge!=0)
            if (xplus_edge in e_sc)
                v1 = e[xplus_edge,1];
                v2 = e[xplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += (G3*invLambda2[xplus_edge] + G2*rho_avg)*(ly*lz/lx)*phixp_past1[xplus_edge];
                divJ_temp += -ecurrent[xplus_edge];
            end
        end
        
        if (yminus_edge != 0)
            if (yminus_edge in e_sc)
                v1 = e[yminus_edge,1];
                v2 = e[yminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += -(G3*invLambda2[yminus_edge] + G2*rho_avg)*(lx*lz/ly)*phiyp_past1[yminus_edge-Ne_x];
                divJ_temp += ecurrent[yminus_edge];
            end
        end
        
        if (yplus_edge != 0)
            if (yplus_edge in e_sc)
                v1 = e[yplus_edge,1];
                v2 = e[yplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += (G3*invLambda2[yplus_edge] + G2*rho_avg)*(lx*lz/ly)*phiyp_past1[yplus_edge-Ne_x];
                divJ_temp += -ecurrent[yplus_edge];
            end
        end

        if (zminus_edge != 0)
            if (zminus_edge in e_sc)
                v1 = e[zminus_edge,1];
                v2 = e[zminus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += -(G3*invLambda2[zminus_edge] + G2*rho_avg)*(lx*ly/lz)*phizp_past1[zminus_edge-Ne_x-Ne_y];
                divJ_temp += ecurrent[zminus_edge];
            end
        end
        
        if (zplus_edge != 0)
            if (zplus_edge in e_sc)
                v1 = e[zplus_edge,1];
                v2 = e[zplus_edge,2];
                rho_avg = 0.5*(rho_past1[v1]+rho_past1[v2]);
                rho_temp  += (G3*invLambda2[zplus_edge] + G2*rho_avg)*(lx*ly/lz)*phizp_past1[zplus_edge-Ne_x-Ne_y];
                divJ_temp += -ecurrent[zplus_edge];
            end
        end
        
        rho[i] = rho_temp*temp + rho_past1[i] + divJ_temp*temp;
    end
    
    return rho;    
end