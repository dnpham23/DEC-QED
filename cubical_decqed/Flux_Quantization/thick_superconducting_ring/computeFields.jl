function computePhi(phix_past1::Array{Float64,1}, phiy_past1::Array{Float64,1}, phixp_current::Array{Float64,1}, phiyp_current::Array{Float64,1}, phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, psixpast1::Array{Float64,1}, psiypast1::Array{Float64,1}, S4::Float64,dt::Float64,ne_x::Int64, Nx::Int64, Ne_x::Int64,Ne_y::Int64, ebound_all::Array{Int64,1}, lx::Float64,ly::Float64)
    phix =  Array{Float64,1}(undef, Ne_x);
    phiy =  Array{Float64,1}(undef, Ne_y);
    for xind=1:Ne_x
        if !(xind in ebound_all)
            rownum   = div(xind-1,ne_x) + 1; # global column
            colnum   = (xind-1)%ne_x + 1;    # global row
            xind_left   = xind - 1;
            y_edge      = (rownum-1)*Nx + colnum; 
            y_edge_left = y_edge-1; 
            phix[xind]  = phix_past1[xind] + (phixp_current[xind]-phixp_past1[xind]) - S4*(dt/2)*( (phixp_past1[xind]^2 - phixp_past1[xind_left]^2)/(lx^2) + (phiyp_past1[y_edge]^2 - phiyp_past1[y_edge_left]^2)/(ly^2) ) - dt*psixpast1[xind];
        else
            phix[xind]   = 0.0;
        end
    end
    
    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        if !(real_yind in ebound_all)
            rownum  = div(yind-1,Nx) + 1;
            colnum  = (yind-1)%Nx + 1;
            y_edge_below = yind - Nx;
            x_edge       = ne_x*(rownum-1) + colnum;
            x_edge_below = x_edge-ne_x; 
            phiy[yind] = phiy_past1[yind] + (phiyp_current[yind]-phiyp_past1[yind]) - S4*(dt/2)*( (phixp_past1[x_edge]^2 - phixp_past1[x_edge_below]^2)/(lx^2) + (phiyp_past1[yind]^2 - phiyp_past1[y_edge_below]^2)/(ly^2) ) - dt*psiypast1[yind];
        else
            phiy[yind]   = 0.0;
        end
    end
    return phix, phiy;
end

function computeBfield3D(phix::Array{Float64,1},phiy::Array{Float64,1},phiz::Array{Float64,1},Ne::Int64, Nx::Int64, Ny::Int64, Nz::Int64,ne_x::Int64, ne_y::Int64, ne_z::Int64, lx::Float64, ly::Float64, lz::Float64,Nex_xyplane::Int64, Ney_xyplane::Int64, Nv_xyplane::Int64)
    N_Bx = ne_y*ne_z*Nx;
    N_By = ne_x*ne_z*Ny;
    N_Bz = ne_x*ne_y*Nz;
    Bx = zeros(Nx, ne_y, ne_z);
    By = zeros(ne_x, Ny, ne_z);
    Bz = zeros(ne_x, ne_y, Nz);
    
    # first, calculate Bx
    for i = 1:Nx
        for jy = 1:ne_y
            for jz = 1:ne_z
                ybelow_edge  = (jz-1)*Ney_xyplane + (jy-1)*Nx + i;
                yabove_edge  = ybelow_edge + Ney_xyplane;
                zbehind_edge = (jz-1)*Nv_xyplane + (jy-1)*Nx + i;
                zfront_edge  = zbehind_edge + Nx;
                Bx[i,jy,jz]   = (phiy[ybelow_edge] - phiy[yabove_edge] + phiz[zfront_edge] - phiz[zbehind_edge])/(ly*lz);
            end
        end
    end
    
    # calculate By
    for i = 1:Ny
        for jx = 1:ne_x
            for jz = 1:ne_z
                xbelow_edge = (jz-1)*Nex_xyplane + (i-1)*ne_x + jx;
                xabove_edge = xbelow_edge + Nex_xyplane;
                zleft_edge  = (jz-1)*Nv_xyplane + (i-1)*Nx + jx;
                zright_edge = zleft_edge + 1;
                By[jx,i,jz]  = (phiz[zleft_edge] - phiz[zright_edge] + phix[xabove_edge] - phix[xbelow_edge])/(lx*lz);
            end
        end
    end
    
    # calculate Bz
    for i = 1:Nz
        for jx = 1:ne_x
            for jy = 1:ne_y
                xbehind_edge = (i-1)*Nex_xyplane + (jy-1)*ne_x + jx;
                xfront_edge  = xbehind_edge + ne_x;
                yleft_edge   = (i-1)*Ney_xyplane + (jy-1)*Nx + jx;
                yright_edge  = yleft_edge + 1;
                Bz[jx,jy,i]  = (phix[xbehind_edge] - phix[xfront_edge] + phiy[yright_edge] - phiy[yleft_edge])/(lx*ly);
            end
        end
    end
    

    return Bx ,By, Bz;
end

function computeCurrent(phixp::Array{Float64,1},phiyp::Array{Float64,1},psix::Array{Float64,1},psiy::Array{Float64,1}, ecurrent::Array{Float64,2}, lx::Float64, ly::Float64, invLambda2::Array{Float64,1}, Is1::Float64, Is2::Float64,ebound_all::Array{Int64,1},Ne_x::Int64, Ne_y::Int64, Nx::Int64, ne_x::Int64, lo::Float64)
    Ix = zeros(Ne_x,1);
    Iy = zeros(Ne_y,1);

    for xind=1:Ne_x
        if !(xind in ebound_all)
            rownum   = div(xind-1,ne_x) + 1; # global column
            colnum   = (xind-1)%ne_x + 1;    # global row
            y_edge_belowleft  = (rownum-2)*Nx + colnum; # index of the yedge is to the left and one row below the current xedge
            y_edge_aboveleft  = y_edge_belowleft+Nx;# the index (in y-edge sc array) of the yedge above and to the left
            Ix[xind] = Is1*(-invLambda2[xind]/mu_list[xind])*phixp[xind]/lx + Is2*eps_list[xind]*( (psix[xind]-psix[xind-1])/lx^2
                                      +(psiy[y_edge_aboveleft]-psiy[y_edge_belowleft])/ly^2 )*phixp[xind]/lx + lo*ecurrent[xind]/ly;
        end
    end
    
    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        if !(real_yind in ebound_all)
            rownum   = div(yind-1,Nx) + 1;
            colnum   = (yind-1)%Nx + 1;
            x_edge_left  = ne_x*(rownum-1) + colnum-1; # the index (in x-edge sc array) of the xedge to the left
            x_edge_right = x_edge_left+1; # the index (in x-edge sc array) of the xedge to the right
            y_edge_below = yind-Nx;
            Iy[yind] = Is1*(-invLambda2[real_yind]/mu_list[real_yind])*phiyp[yind]/ly + Is2*eps_list[real_yind]*( 
                                (psix[x_edge_right] - psix[x_edge_left])/lx^2 +(psiy[yind]-psiy[y_edge_below])/ly^2 )*phiyp[yind]/ly + lo*ecurrent[real_yind]/lx;
        end
    end
    return Ix, Iy;
end


function computeCharge(psix::Array{Float64,1},psiy::Array{Float64,1}, ecurrent::Array{Float64,2}, lx::Float64, ly::Float64, eps::Float64, mu::Float64, lambda::Float64, Nx::Int64, Nv::Int64, rho1::Float64, v2exmap::Array{Int,2}, v2eymap::Array{Int,2},v_sc::Array{Int64,1}, Ne_x::Int64)
    
    Nv_sc_all = length(v_sc);
    rho       = zeros(Nv,1);
    for i=1:Nv_sc_all
        ind = v_sc[i];
        xminus_edge = v2exmap[ind,1];
        xplus_edge  = v2exmap[ind,2];
        yminus_edge = v2eymap[ind,1];
        yplus_edge  = v2eymap[ind,2];
        temp = 0.0;
        if (xminus_edge!=0)
            temp += -psix[xminus_edge]/lx^2;
        end
        
        if (xplus_edge!=0)
            temp += psix[xplus_edge]/lx^2;
        end
        
        if (yminus_edge!=0)
            temp += -psiy[yminus_edge-Ne_x]/ly^2;
        end
        
        if (yplus_edge!=0)
            temp += psiy[yplus_edge-Ne_x]/ly^2;
        end
        rho[ind] = -temp*rho1*(eps*mu*lambda^2);
    end
    return rho;
end

function computeCurrent_JJ3D(phixp::Array{Float64,1},phiyp::Array{Float64,1}, phizp::Array{Float64,1}, rho_list::Array{Float64,1}, ecurrent::Array{Float64,2}, ecurrentJJ::Array{Float64,2}, lx::Float64, ly::Float64, lz::Float64, invLambda2::Array{Float64,1}, Is1::Float64, Is3::Float64,ebound_all::Array{Int64,1}, Ne_x::Int64, Ne_y::Int64, Ne_z::Int64, Nx::Int64, ne_x::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64, Nv_xyplane::Int64)
    Ix = zeros(Ne_x,1);
    Iy = zeros(Ne_y,1);
    Iz = zeros(Ne_z,1);

    for xind=1:Ne_x
        if !(xind in ebound_all)
            ie_xyplane = (xind-1)%Nex_xyplane + 1;
            plane_ind  = div(xind-1,Nex_xyplane) + 1;
            rownum     = div(ie_xyplane-1,ne_x) + 1; # column in the xy plane
            colnum     = (ie_xyplane-1)%ne_x + 1;    # row in the xy plane
            vind_left  = (plane_ind-1)*Nv_xyplane + (rownum-1)*Nx + colnum;
            vind_right = vind_left + 1;
            rho_avg    = 0.5*(rho_list[vind_left] + rho_list[vind_right]);
            
            #Ix[xind] = -(Is1*invLambda2[xind] + Is3*rho_avg)*phixp[xind]/lx + ecurrent[xind] + ecurrentJJ[xind];
            Ix[xind] = -(Is1*invLambda2[xind] + Is3*rho_avg)*phixp[xind]/lx + ecurrentJJ[xind];
        end
    end
    
    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        if !(real_yind in ebound_all)
            ie_xyplane  = (yind-1)%Ney_xyplane + 1;
            plane_ind   = div(yind-1,Ney_xyplane) + 1;
            rownum      = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
            colnum      = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
            vind_ybelow = (plane_ind-1)*Nv_xyplane + (rownum-1)*Nx + colnum;
            vind_yabove = vind_ybelow + Nx;
            rho_avg     = 0.5*(rho_list[vind_ybelow] + rho_list[vind_yabove]);
      
            #Iy[yind] = -(Is1*invLambda2[real_yind] + Is3*rho_avg)*phiyp[yind]/ly + ecurrent[real_yind] + ecurrentJJ[real_yind];
            Iy[yind] = -(Is1*invLambda2[real_yind] + Is3*rho_avg)*phiyp[yind]/ly + ecurrentJJ[real_yind];
        end
    end
    
    for zind=1:Ne_z
        real_zind = zind + Ne_x + Ne_y;
        if !(real_zind in ebound_all)
            ie_xyplane  = (zind-1)%Nv_xyplane + 1;
            plane_ind   = div(zind-1,Nv_xyplane) + 1;
            rownum      = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
            colnum      = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
            vind_zbelow = zind;
            vind_zabove = vind_zbelow + Nv_xyplane;
            rho_avg     = 0.5*(rho_list[vind_zbelow] + rho_list[vind_zabove]);
            #Iz[zind] = -(Is1*invLambda2[real_zind] + Is3*rho_avg)*phizp[zind]/lz + ecurrent[real_zind] + ecurrentJJ[real_zind];
            Iz[zind] = -(Is1*invLambda2[real_zind] + Is3*rho_avg)*phizp[zind]/lz + ecurrentJJ[real_zind];
        end
    end
    return Ix, Iy, Iz;
end
