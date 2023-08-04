function reducedMat(tree::Array{Int,1},Ne::Int64,Nbranch::Int64,ne_x::Int64,ne_y::Int64,Ne_x::Int64,Nx::Int64)
    #x = Array{Float64,2}(undef, Ne, nbranch);
    x = zeros(Ne, Nbranch);
    lastbranch = tree[Nbranch];
    for i = 1:Ne
        isbranch = findfirst(isequal(i),tree);
        colnum = (i-1)%ne_x + 1;
        if  isbranch != nothing  # if this edge is a branch of tree
            x[i,isbranch] = 1;
        elseif lastbranch < i # if this edge is in the last row of branches (that was not in tree)
            colnum = (i-Ne_x-1)%Nx + 1; 
            x[i,colnum] = 1; # this is to impose zero BC for this edge
            #for j = 1:(ne_y-1)
            #    iremain = (i-Ne_x-1)%Nx + 1;
            #    btemp   = ne_x + (j-1)*Nx + iremain;
            #    x[i,btemp] = -1;
            #end
        elseif (i < (Ne_x-ne_x))&&(colnum != 1)&&(colnum != ne_x)   # if this edge is not a branch, and it is not at the boundary
            rownum = div(i-1,ne_x) + 1;
            colnum = (i-1)%ne_x + 1;
            for j = 1:(rownum-1)
                btemp = ne_x + (j-1)*Nx + colnum;
                x[i,btemp]   = -1;
                x[i,btemp+1] = 1;
            end
        else # the remaining scenerio is this edge is at the boundary, and it's not a branch
            colnum = (i-1)%ne_x + 1;
            x[i,colnum] = 1;
        end
    end
    return x;
end


function GaussMat(v2exmap::Array{Int,2},v2eymap::Array{Int,2},Nv::Int64,Ne::Int64,Nx::Int64,lx::Float64,ly::Float64)
    Y = zeros(Nv,Ne);
    xweight = ly/lx;
    yweight = lx/ly;
    for i = 1:Nv
        xminus_edge = v2exmap[i,1];
        xplus_edge  = v2exmap[i,2];
        yminus_edge = v2eymap[i,1];
        yplus_edge  = v2eymap[i,2];
        if xminus_edge!=0
            Y[i,xminus_edge] = -xweight;
        end
        if xplus_edge!=0
            Y[i,xplus_edge]  = xweight;
        end
        if yminus_edge!=0
            Y[i,yminus_edge] = -yweight;
        end
        if yplus_edge!=0
            Y[i,yplus_edge]  = yweight;
        end
    end
    vredund = [div(Nx,2)+1 [Nv-Nx+1:Nv;]'];
    Y = Y[setdiff(1:end, vredund), :];
    return Y;
end

function charge(v::Array{Float64,2},chargeloc::Array{Float64,2},chargeamp::Array{Float64,1},xtol::Float64,ytol::Float64,Nv::Int64,Nx::Int64, eps_air::Float64,eps_sc::Float64,psixpast1::Array{Float64,1},psiypast1::Array{Float64,1},phixp_past1::Array{Float64,1},phiyp_past1::Array{Float64,1}, v2exmap::Array{Int,2},v2eymap::Array{Int,2},v_sc::Array{Int,1},Nv_sc::Int64,A2::Float64,G1::Float64,G2::Float64,dt::Float64, mu_sc::Float64, lambda::Float64,lx::Float64,ly::Float64, e_scx::Array{Int,1}, e_scy::Array{Int,1},ne_scx_seg::Int64, ne_scy_seg::Int64, v_interface::Array{Int,1},Ne_x::Int64)
    #ytol = xtol*1.0;
    numcharge = length(chargeloc[:,1]);
    vcharge = zeros(Nv,1);
    for i = 1:numcharge
        for j = 1:Nv
            if (abs(chargeloc[i,1]-v[j,1])<xtol)&&(abs(chargeloc[i,2]-v[j,2])<ytol)
                vcharge[j] = A2*chargeamp[i]/eps_air;
                break;
            end
        end
    end

    vredund = [div(Nx,2)+1 [Nv-Nx+1:Nv;]'];
    vcharge = vcharge[setdiff(1:end, vredund), :];
    return vcharge;
end


function charge2(v::Array{Float64,2},chargeloc::Array{Float64,2},chargeamp::Array{Float64,1},xtol::Float64,ytol::Float64,Nv::Int64, Nx::Int64, eps_air::Float64,eps_sc::Float64,psixpast1::Array{Float64,1},psiypast1::Array{Float64,1},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, v2exmap::Array{Int,2},v2eymap::Array{Int,2},A2::Float64,G1::Float64, G2::Float64, dt::Float64, mu_sc::Float64, lambda::Float64,lx::Float64,ly::Float64, Ne_x::Int64,v_sc::Array{Int,1})
    #ytol = xtol*1.0;
    numcharge = length(chargeloc[:,1]);
    vcharge = zeros(Nv,1);
    for i = 1:numcharge
        for j = 1:Nv
            if (abs(chargeloc[i,1]-v[j,1])<xtol)&&(abs(chargeloc[i,2]-v[j,2])<ytol)
                vcharge[j] = A2*chargeamp[i]/eps_air;
                break;
            end
        end
    end

    # Compute the 'charges' at the vertices inside the superconductor
    # These charges are actually functions of the fields and field gradients surrounding the vertex
    # There are many ways to do this, but first let's consider a simple approximation for the gradient
    for i = 1:Nv
        xminus_edge = v2exmap[i,1];
        xplus_edge  = v2exmap[i,2];
        yminus_edge = v2eymap[i,1];
        yplus_edge  = v2eymap[i,2];

        temp1 = 0.0;
        temp2 = 0.0;
        
        if (i in v_sc)
            w = -G1/(eps_sc*mu_sc*lambda^2);
        else
            w = 0.0;
        end
        if (xminus_edge!=0)
            w     = w - G2*psixpast1[xminus_edge]/lx^2;
            temp1 = temp1 - (ly/lx)*psixpast1[xminus_edge]; # psi is indexed globally
            temp2 = temp2 - (ly/lx)*phixp_past1[xminus_edge]; # phi is index in the sc region
        end

        if (xplus_edge!=0)
            w     = w + G2*psixpast1[xplus_edge]/lx^2;
            temp1 = temp1 + (ly/lx)*psixpast1[xplus_edge]; # psi is indexed globally
            temp2 = temp2 + (ly/lx)*phixp_past1[xplus_edge]; # phi is index in the sc region
        end

        if (yminus_edge!=0)
            w     = w - G2*psiypast1[yminus_edge-Ne_x]/ly^2;
            temp1 = temp1 - (lx/ly)*psiypast1[yminus_edge-Ne_x]; # psi is indexed globally
            temp2 = temp2 - (lx/ly)*phiyp_past1[yminus_edge-Ne_x]; # phi is index in the sc region
        end

        if (yplus_edge!=0)
            w     = w + G2*psiypast1[yplus_edge-Ne_x]/ly^2;
            temp1 = temp1 + (lx/ly)*psiypast1[yplus_edge-Ne_x]; # psi is indexed globally
            temp2 = temp2 + (lx/ly)*phiyp_past1[yplus_edge-Ne_x]; # phi is index in the sc region
        end
        
        if vcharge[i] == 0.0
            vcharge[i] += dt*w*temp2 + temp1;
        end
    end

    vredund = [div(Nx,2)+1 [Nv-Nx+1:Nv;]'];
    vcharge = vcharge[setdiff(1:end, vredund), :];
    return vcharge;
end


function charge3(v::Array{Float64,2},chargeloc::Array{Float64,2},chargeamp::Array{Float64,1},xtol::Float64,ytol::Float64,Nv::Int64, Nx::Int64, eps_air::Float64,eps_sc::Float64,psixpast1::Array{Float64,1},psiypast1::Array{Float64,1},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, v2exmap::Array{Int,2},v2eymap::Array{Int,2},v_sc::Array{Int,1},Nv_sc::Int64,A2::Float64,G1::Float64, G2::Float64, dt::Float64, mu_sc::Float64, lambda::Float64,lx::Float64,ly::Float64, e_scx::Array{Int,1}, e_scy::Array{Int,1},ne_scx_seg::Int64, ne_scy_seg::Int64, v_interface::Array{Int,1},Ne_x::Int64)
    #ytol = xtol*1.0;
    numcharge = length(chargeloc[:,1]);
    vcharge = zeros(Nv,1);
    for i = 1:numcharge
        for j = 1:Nv
            if (abs(chargeloc[i,1]-v[j,1])<xtol)&&(abs(chargeloc[i,2]-v[j,2])<ytol)
                vcharge[j] = A2*chargeamp[i]/eps_air;
                break;
            end
        end
    end

    # Compute the 'charges' at the vertices inside the superconductor
    # These charges are actually functions of the fields and field gradients surrounding the vertex
    # There are many ways to do this, but first let's consider a simple approximation for the gradient
    for i = 1:Nv
        xminus_edge = v2exmap[i,1];
        xplus_edge  = v2exmap[i,2];
        yminus_edge = v2eymap[i,1];
        yplus_edge  = v2eymap[i,2];

        temp1 = 0.0;
        temp2 = 0.0;
        
        if (i in v_sc)
            w = -G1/(eps_sc*mu_sc*lambda^2);
        else
            w = 0.0;
        end
        if (xminus_edge!=0)&&(xplus_edge!=0)
            w     = w - G2*psixpast1[xminus_edge]/lx^2;
            temp1 = temp1 - (ly/lx)*psixpast1[xminus_edge]; # psi is indexed globally
            temp2 = temp2 - (ly/lx)*phixp_past1[xminus_edge]; # phi is index in the sc region

            w     = w + G2*psixpast1[xplus_edge]/lx^2;
            temp1 = temp1 + (ly/lx)*psixpast1[xplus_edge]; # psi is indexed globally
            temp2 = temp2 + (ly/lx)*phixp_past1[xplus_edge]; # phi is index in the sc region
        end

        if (yminus_edge!=0)&&(yplus_edge!=0)
            w     = w - G2*psiypast1[yminus_edge-Ne_x]/ly^2;
            temp1 = temp1 - (lx/ly)*psiypast1[yminus_edge-Ne_x]; # psi is indexed globally
            temp2 = temp2 - (lx/ly)*phiyp_past1[yminus_edge-Ne_x]; # phi is index in the sc region

            w     = w + G2*psiypast1[yplus_edge-Ne_x]/ly^2;
            temp1 = temp1 + (lx/ly)*psiypast1[yplus_edge-Ne_x]; # psi is indexed globally
            temp2 = temp2 + (lx/ly)*phiyp_past1[yplus_edge-Ne_x]; # phi is index in the sc region
        end
        
        if vcharge[i] == 0.0
            vcharge[i] += dt*w*temp2 + temp1;
        end
    end

    vredund = [div(Nx,2)+1 [Nv-Nx+1:Nv;]'];
    vcharge = vcharge[setdiff(1:end, vredund), :];
    return vcharge;
end


function current3D(v::Array{Float64,2},currentloc::Array{Float64,3},currentamp::Array{Float64,1},xtol::Float64,ytol::Float64,ztol::Float64, Nv::Int64, Nx::Int64, Ne::Int64,e::Array{Int,2})
    numcurrent = length(currentamp);
    ecurrent = zeros(Ne,1);
    v1 = 0;
    v2 = 0;
    for i = 1:numcurrent
        x1 = currentloc[1,1,i];
        y1 = currentloc[1,2,i];
        z1 = currentloc[1,3,i];
        x2 = currentloc[2,1,i];
        y2 = currentloc[2,2,i];
        z2 = currentloc[2,3,i];
        v1 = 0;
        v2 = 0;
        for j = 1:Nv
            if (abs(x1-v[j,1])<xtol)&&(abs(y1-v[j,2])<ytol)&&(abs(z1-v[j,3])<ztol)
                v1 = j;
            end
            if (abs(x2-v[j,1])<xtol)&&(abs(y2-v[j,2])<ytol)&&(abs(z2-v[j,3])<ztol)
                v2 = j;
            end
            if (v1 != 0)&&(v2 != 0)
                break;
            end
        end
        for k = 1:Ne
            if ((v1==e[k,1])&&(v2==e[k,2]))||(v2==e[k,1])&&(v1==e[k,2])
                ecurrent[k] = currentamp[i];
                break;
            end
        end
    end
    return ecurrent;
end

function Ampere_AMat(S4::Float64,dt::Float64,Ne::Int64,ne_x::Int64, Nx::Int64,Ne_x::Int64,Ne_y::Int64,ebound::Array{Int64,1}, e_sc::Array{Int,1}, lx::Float64,ly::Float64,Ny::Int64,ne_y::Int64, eps_list::Array{Float64,1}, mu_list::Array{Float64,1})

    AMat = zeros(Ne, Ne);
    # x-edges
    for xind = 1:Ne_x
        if !(xind in ebound)
            rownum  = div(xind-1,ne_x) + 1; # global row
            colnum  = (xind-1)%ne_x + 1;    # global column
            yind    = Ne_x + (rownum-1)*Nx + colnum;    # actual index in the global grid
            if (colnum != 1)&&(colnum != ne_x)
                # if this x edge doesn't touch the vertical boundary on the left
                AMat[xind,xind]   = -S4*(mu_list[xind]*eps_list[xind])/(2.0*dt)/lx^2;
                AMat[xind,yind]   = -S4*(mu_list[yind]*eps_list[yind])/(2.0*dt)/ly^2;
                AMat[xind,xind-1] =  S4*(mu_list[xind-1]*eps_list[xind-1])/(2.0*dt)/lx^2;
                AMat[xind,yind-1] =  S4*(mu_list[yind-1]*eps_list[yind-1])/ly^2;
            end
        end
    end
    # y-edges
    for i = 1:Ne_y
        yind = i + Ne_x;
        if !(yind in ebound)
            rownum = div(yind-Ne_x-1,Nx) + 1;
            colnum = (yind-Ne_x-1)%Nx + 1;
            xind   = (rownum-1)*ne_x + colnum;
            if (rownum != ne_y)&&(rownum != 1)
                # if the bottom of this y edge doesn't touch a horizontal boundary
                AMat[yind,yind]      = -S4*(mu_list[yind]*eps_list[yind])/(2.0*dt)/ly^2;
                AMat[yind,xind]      = -S4*(mu_list[xind]*eps_list[xind])/(2.0*dt)/lx^2;
                AMat[yind,yind-Nx]   =  S4*(mu_list[yind-Nx]*eps_list[yind-Nx])/(2.0*dt)/ly^2;
                AMat[yind,xind-ne_x] =  S4*(mu_list[xind-ne_x]*eps_list[xind-ne_x])/(2.0*dt)/lx^2;
            end
        end
    end
    return AMat;
end

function Ampere_BMat(dt::Float64,Ne::Int64,Nx::Int64, Ny::Int64,ne_x::Int64,ebound_all::Array{Int64,1},eps_list::Array{Float64,1}, mu_list::Array{Float64,1})
    
    BMat = zeros(Ne, Ne);    
    for i = 1:Ne
        # only fill in the diagonal element if the edge is not at boundary
        if !(i in ebound_all)
            BMat[i,i] =  mu_list[i]*eps_list[i]/(dt^2);
        end
    end
    return BMat;
end

function Ampere_Ccol(J1::Float64,S1::Float64,S2::Float64,S3::Float64,S4::Float64,dt::Float64,ne_x::Int64,Nx::Int64,Ne_x::Int64,Ne_y::Int64, ebound_all::Array{Int64,1},psixpast1::Array{Float64,1},psiypast1::Array{Float64,1},phixp_past1::Array{Float64,1},phiyp_past1::Array{Float64,1}, phixp_past2::Array{Float64,1},phiyp_past2::Array{Float64,1}, e_sc::Array{Int,1},Ny::Int64,e_scx::Array{Int,1},e_scy::Array{Int,1}, lx::Float64,ly::Float64,ne_y::Int64,eps_list::Array{Float64,1}, mu_list::Array{Float64,1}, invLambda2::Array{Float64,1},Ne::Int64, ecurrent::Array{Float64,2})

    Ccol = zeros(Ne,1);

    for xind=1:Ne_x
        if !(xind in ebound_all)
            rownum   = div(xind-1,ne_x) + 1; # global column
            colnum   = (xind-1)%ne_x + 1;    # global row
            xedge_below = xind-ne_x; # find the index (in x-edge sc array) of the xedge below current xedge
            xedge_above = xind+ne_x; # find the index (in x-edge sc array) of the xedge above current xedge
            y_edge_belowleft  = (rownum-2)*Nx + colnum; # index of the yedge is to the left and one row below the current xedge
            y_edge_belowright = y_edge_belowleft+1; # the index (in y-edge sc array) of the yedge below and to the right
            y_edge_aboveleft  = y_edge_belowleft+Nx;# the index (in y-edge sc array) of the yedge above and to the left
            y_edge_aboveright = y_edge_belowleft+Nx+1;# the index (in y-edge sc array) of the yedge above and to the right

            temp1    = (S1/ly^2)*(-phixp_past1[xedge_below] + 2*phixp_past1[xind] -phixp_past1[xedge_above] + phiyp_past1[y_edge_belowleft] - phiyp_past1[y_edge_belowright] - phiyp_past1[y_edge_aboveleft] + phiyp_past1[y_edge_aboveright]);
            temp2    = S2*invLambda2[xind]*phixp_past1[xind];
            leftover = -mu_list[xind]*eps_list[xind]/(dt^2)*(2*phixp_past1[xind]-phixp_past2[xind]);
            xedge_left        = xind-1; # find the index (in x-edge sc array) of the xedge to the left
            y_edge_aboveleft2 = y_edge_belowleft+Nx-1;# index (in yedge sc array) of yedge above and twice to the left
            temp3 = -S3*eps_list[xind]*( (psixpast1[xind]-psixpast1[xind-1])/lx^2
                                      +(psiypast1[y_edge_aboveleft]-psiypast1[y_edge_belowleft])/ly^2 )*phixp_past1[xind];
            temp4 = S4*(mu_list[xind]*eps_list[xind])/(2.0*dt)*( (phixp_past1[xind]^2 - phixp_past1[xedge_left]^2)/lx^2
                                                            +(phiyp_past1[y_edge_aboveleft]^2 - phiyp_past1[y_edge_aboveleft2]^2)/ly^2 );
            Ccol[xind] = temp1+temp2+temp3+temp4+leftover - J1*mu_list[xind]*lx/ly*ecurrent[xind];                
        end
    end


    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        if !(real_yind in ebound_all)
            rownum   = div(yind-1,Nx) + 1;
            colnum   = (yind-1)%Nx + 1;
            y_edge_left  = yind-1; # the index (in y-edge sc array) of the yedge to the left
            y_edge_right = yind+1; # the index (in y-edge sc array) of the yedge to the right
            x_edge_left  = ne_x*(rownum-1) + colnum-1; # the index (in x-edge sc array) of the xedge to the left
            x_edge_right = x_edge_left+1; # the index (in x-edge sc array) of the xedge to the right
            x_edge_left_above  = x_edge_left+ne_x; # the index (in x-edge sc array) of the xedge above and to the left
            x_edge_right_above = x_edge_left+ne_x+1; # the index (in x-edge sc array) of the xedge above to the right

            temp1    = (S1/lx^2)*(-phiyp_past1[y_edge_left] + 2*phiyp_past1[yind] - phiyp_past1[y_edge_right] + phixp_past1[x_edge_left] - phixp_past1[x_edge_left_above] - phixp_past1[x_edge_right] + phixp_past1[x_edge_right_above]);
            temp2    = S2*invLambda2[real_yind]*phiyp_past1[yind];
            leftover = -mu_list[real_yind]*eps_list[real_yind]/(dt^2)*(2*phiyp_past1[yind]-phiyp_past2[yind]);

            # if the bottom of this y edge doesn't touch the horizontal boundary
            y_edge_below       = yind-Nx; # the index (in y-edge sc array) of the yedge below current yedge
            x_edge_right_below = x_edge_left-ne_x+1; # index (in y-edge sc array) of xedge below to the right
            temp3 = -S3*eps_list[real_yind]*( (psixpast1[x_edge_right]-psixpast1[x_edge_left])/lx^2
                                       +(psiypast1[yind]-psiypast1[y_edge_below])/ly^2 )*phiyp_past1[yind];
            temp4 = S4*(eps_list[real_yind]*mu_list[real_yind])/(2.0*dt)*( (phixp_past1[x_edge_right]^2 
                    - phixp_past1[x_edge_right_below]^2 )/lx^2 +(phiyp_past1[yind]^2 - phiyp_past1[y_edge_below]^2)/ly^2 );
            Ccol[real_yind] = temp1+temp2+temp3+temp4+leftover - J1*mu_list[real_yind]*ly/lx*ecurrent[real_yind];
        end
    end
    
    return Ccol;
end

using LinearAlgebra

function Ampere_MMat_floatingBC(AMat::Array{Float64,2},BMat::Array{Float64,2},Ccol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1},ebound_all::Array{Int64,1},Ne::Int64, Ne_x::Int64,ne_x::Int64,Nx::Int64,Ny::Int64)
    MMat = zeros(Ne,Ne);
    phip_past1      = [phixp_past1; phiyp_past1];
    phip_past1_diag = Diagonal(vec(phip_past1));
    MMat = 2.0*AMat*phip_past1_diag + BMat;
    
    for i=1:length(ebound_all)
        ind = ebound_all[i];
        if ind <= Ne_x
            rownum   = div(ind-1,ne_x) + 1; # global column
            colnum   = (ind-1)%ne_x + 1;    # global row
            if (rownum!= 1)&&(rownum != Ny)
                if (colnum==1)
                    MMat[ind,ind] = 1.0;
                    MMat[ind,ind+1] = -1.0;
                else
                    MMat[ind,ind] = 1.0;
                    MMat[ind, ind-1] = -1.0;
                end
            elseif (rownum == 1)
                MMat[ind,ind] = 1.0;
                MMat[ind,ind+ne_x] = -1.0;
            else
                MMat[ind,ind] = 1.0;
                MMat[ind,ind-ne_x] = -1.0;
            end
        else
            rownum = div(ind-1-Ne_x,Nx) + 1;
            colnum = (ind-1-Ne_x)%Nx + 1;
            if (colnum != 1)&&(colnum != Nx)
                if (rownum==1)
                    MMat[ind,ind] = 1.0;
                    MMat[ind,ind+Nx] = -1.0;
                else
                    MMat[ind,ind] = 1.0;
                    MMat[ind,ind-Nx] = -1.0;
                end
            elseif (colnum == 1)
                MMat[ind,ind] = 1.0;
                MMat[ind,ind+1] = -1.0;
            else
                MMat[ind,ind] = 1.0;
                MMat[ind,ind-1] = -1.0;
            end
        end
    end
    return MMat;
end

function Ampere_MMat(AMat::Array{Float64,2},BMat::Array{Float64,2},Ccol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1},ebound_all::Array{Int64,1}, Ne::Int64)
    MMat = zeros(Ne,Ne);
    phip_past1      = [phixp_past1; phiyp_past1];
    phip_past1_diag = Diagonal(vec(phip_past1));
    MMat = 2.0*AMat*phip_past1_diag + BMat;
    
    for i=1:length(ebound_all)
        ind = ebound_all[i];
        MMat[ind,ind] = 1.0;
    end

    return MMat;
end

function Ampere_Kcol(AMat::Array{Float64,2},BMat::Array{Float64,2},Ccol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1})
    phip_past1      = [phixp_past1; phiyp_past1];
    Kcol = AMat*(phip_past1.*phip_past1) + BMat*phip_past1 + Ccol;
    return Kcol;
end

function SteadyState_HMat_3D(Is1::Float64,Ne::Int64,ne_x::Int64,Nx::Int64,Ne_x::Int64,Ne_y::Int64,ebound::Array{Int64,1}, lx::Float64,ly::Float64,lz::Float64,Ny::Int64,ne_y::Int64,invLambda2::Array{Float64,1},Nex_xyplane::Int64, Ney_xyplane::Int64, Nv_xyplane::Int64)
    HMat = zeros(Ne,Ne);
    for xind=1:Ne_x
        if !(xind in ebound)
            ie_xyplane = (xind-1)%Nex_xyplane + 1;
            plane_ind  = div(xind-1,Nex_xyplane) + 1;
            rownum     = div(ie_xyplane-1,ne_x) + 1; # column in the xy plane
            colnum     = (ie_xyplane-1)%ne_x + 1;    # row in the xy plane
            
            xedge_y_below     = xind - ne_x; # the xedge below (in y direction) the current xedge
            xedge_y_above     = xind + ne_x; # the xedge above (in y direction) the current xedge           
            y_edge_belowleft  = (rownum-2)*Nx + colnum + (plane_ind-1)*Ney_xyplane; # yedge to the left and one row below the current xedge
            y_edge_belowright = y_edge_belowleft + 1; # the yedge below and to the right
            y_edge_aboveleft  = y_edge_belowleft + Nx;# the yedge above and to the left
            y_edge_aboveright = y_edge_belowleft + Nx+1;# the yedge above and to the right
            
            xedge_z_below     = xind - Nex_xyplane; # the xedge below (in z direction) the current xedge
            xedge_z_above     = xind + Nex_xyplane; # the xedge above (in z direction) the current xedge
            z_edge_belowleft  = (rownum-1)*Nx + colnum + (plane_ind-2)*Nv_xyplane; # zedge to the left and one row below the current xedge
            z_edge_belowright = z_edge_belowleft + 1;
            z_edge_aboveleft  = z_edge_belowleft + Nv_xyplane;
            z_edge_aboveright = z_edge_aboveleft + 1;
            
            tempy = 1/ly^2;
            HMat[xind, xedge_y_below] = -tempy;
            HMat[xind, xind]          = 2*tempy + invLambda2[xind];
            HMat[xind, xedge_y_above] = -tempy;
            HMat[xind, y_edge_belowleft+Ne_x]  =  tempy;
            HMat[xind, y_edge_belowright+Ne_x] = -tempy;
            HMat[xind, y_edge_aboveleft+Ne_x]  = -tempy;
            HMat[xind, y_edge_aboveright+Ne_x] =  tempy;
            
            tempz = 1/lz^2;
            HMat[xind, xedge_z_below]  = -tempz;
            HMat[xind, xind]          += 2*tempz;
            HMat[xind, xedge_z_above]  = -tempz;
            HMat[xind, z_edge_belowleft+Ne_x+Ne_y]  =  tempz;
            HMat[xind, z_edge_belowright+Ne_x+Ne_y] = -tempz;
            HMat[xind, z_edge_aboveleft+Ne_x+Ne_y]  = -tempz;
            HMat[xind, z_edge_aboveright+Ne_x+Ne_y] =  tempz;
        end
    end
    
    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        if !(real_yind in ebound)
            ie_xyplane = (yind-1)%Ney_xyplane + 1;
            plane_ind  = div(yind-1,Ney_xyplane) + 1;
            rownum     = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
            colnum     = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
            
            y_edge_left  = yind-1; # the index (in y-edge sc array) of the yedge to the left
            y_edge_right = yind+1; # the index (in y-edge sc array) of the yedge to the right
            x_edge_left  = ne_x*(rownum-1) + colnum-1 + (plane_ind-1)*Nex_xyplane; # xedge to the left
            x_edge_right = x_edge_left + 1; # xedge to the right
            x_edge_left_above  = x_edge_left+ne_x; # the xedge above (in y direction) and to the left
            x_edge_right_above = x_edge_left+ne_x+1; # the xedge above (in y direction) to the right

            y_edge_below = yind - Ney_xyplane;  # the yedge below (in z direction) the current yedge
            y_edge_above = yind + Ney_xyplane;  # the yedge below (in z direction) the current yedge
            z_edge_behind = ie_xyplane + (plane_ind-2)*Nv_xyplane; # zedge behind (in y direction) and one row below the current yedge
            z_edge_front  = z_edge_behind + Nx; # zedge in front (in y direction) and one row below the current yedge
            z_edge_behindabove = z_edge_behind + Nv_xyplane; # zedge behind (in y direction) and above the current yedge
            z_edge_frontabove  = z_edge_behindabove + Nx;    # zedge in front (in y direction) and above the current yedge
            
            tempx = 1/lx^2;
            HMat[real_yind, y_edge_left+Ne_x]   = -tempx;
            HMat[real_yind, real_yind]          =  2*tempx + invLambda2[real_yind];
            HMat[real_yind, y_edge_right+Ne_x]  = -tempx; 
            HMat[real_yind, x_edge_left]        =  tempx;
            HMat[real_yind, x_edge_left_above]  = -tempx;
            HMat[real_yind, x_edge_right]       = -tempx;
            HMat[real_yind, x_edge_right_above] =  tempx;
            
            tempz = 1/lz^2;
            HMat[real_yind, y_edge_below+Ne_x]       = -tempz;
            HMat[real_yind, real_yind]              +=  2*tempz;
            HMat[real_yind, y_edge_above+Ne_x]       = -tempz;
            HMat[real_yind, z_edge_behind+Ne_x+Ne_y] =  tempz;
            HMat[real_yind, z_edge_front+Ne_x+Ne_y]  = -tempz;
            HMat[real_yind, z_edge_behindabove+Ne_x+Ne_y] = -tempz;
            HMat[real_yind, z_edge_frontabove+Ne_x+Ne_y]  =  tempz;
        end
    end
    
    for zind=1:Ne_z
        real_zind = zind + Ne_x + Ne_y;
        if !(real_zind in ebound)
            ie_xyplane = (zind-1)%Nv_xyplane + 1;
            plane_ind   = div(zind-1,Nv_xyplane) + 1;
            rownum     = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
            colnum     = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
            
            z_edge_left  = zind-1; # the zedge to the left
            z_edge_right = zind+1; # the zedge to the right
            x_edge_left  = ne_x*(rownum-1) + colnum-1 + (plane_ind-1)*Nex_xyplane; # xedge to the left
            x_edge_right = x_edge_left + 1; # xedge to the right
            x_edge_left_above  = x_edge_left+Nex_xyplane; # the xedge above (in z direction) and to the left
            x_edge_right_above = x_edge_left_above+1; # the xedge above (in z direction) to the right
            
            z_edge_behind = zind-Nx; # zedge behind (in y direction) the current zedge
            z_edge_front  = zind+Nx; # zedge in front (in y direction) the current zedge
            y_edge_behind = (rownum-2)*Nx + colnum + (plane_ind-1)*Ney_xyplane; # yedge behind (in y direction) the current zedge
            y_edge_front  = y_edge_behind + Nx; # yedge in front (in y direction) the current zedge
            y_edge_behindabove = y_edge_behind + Ney_xyplane; # yedge behind (in y direction) and above the current zedge
            y_edge_frontabove  = y_edge_behindabove + Nx; # yedge in front (in y direction) and above the current zedge
            
            tempx = 1/lx^2;
            HMat[real_zind, z_edge_left+Ne_x+Ne_y]  = -tempx;
            HMat[real_zind, real_zind]              =  2*tempx + invLambda2[real_zind];
            HMat[real_zind, z_edge_right+Ne_x+Ne_y] = -tempx;
            HMat[real_zind, x_edge_left]            =  tempx;
            HMat[real_zind, x_edge_right]           = -tempx;
            HMat[real_zind, x_edge_left_above]      = -tempx;
            HMat[real_zind, x_edge_right_above]     =  tempx;
            
            tempy = 1/ly^2;
            HMat[real_zind, z_edge_behind+Ne_x+Ne_y] = -tempy;
            HMat[real_zind, real_zind]              +=  2*tempy;
            HMat[real_zind, z_edge_front+Ne_x+Ne_y]  = -tempy;
            HMat[real_zind, y_edge_behind+Ne_y]      =  tempy;
            HMat[real_zind, y_edge_front+Ne_y]       = -tempy;
            HMat[real_zind, y_edge_behindabove+Ne_y] = -tempy;
            HMat[real_zind, y_edge_frontabove+Ne_y]  =  tempy;
        end
    end
    for i = 1:length(ebound)
        ind = ebound[i];
        HMat[ind,:] .= 0.0;
        HMat[ind,ind] = 1.0;
    end
    HMat = Is1*HMat;
    return HMat;
end

function SteadyState_Kcol(ecurrent::Array{Float64,2},lx::Float64,ly::Float64, lz::Float64)
    Kcol = zeros(Ne);
    Kcol[1:Ne_x] = lx.*ecurrent[1:Ne_x];
    Kcol[Ne_x+1:Ne_x+Ne_y] = ly.*ecurrent[Ne_x+1:Ne_x+Ne_y];
    Kcol[Ne_x+Ne_y+1:end]  = lz.*ecurrent[Ne_x+Ne_y+1:end];
    return Kcol;
end

function reducedPhiMat3D(tree::Array{Int,1},Ne::Int64,Nbranch::Int64,ne_x::Int64,ne_y::Int64,ne_z::Int64, Ne_x::Int64, Ne_y::Int64, Nx::Int64,Ny::Int64, Nz::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64, Nv_xyplane::Int64)
    Nx_reducedPhi = ne_x*ne_y + ne_x*Ny*ne_z;
    Ny_reducedPhi = ne_y*Nx*ne_z;
    Nreduced = Nx_reducedPhi + Ny_reducedPhi;
    reduced_ind_map = zeros(Nreduced);
    z = zeros(Ne,Nreduced);
    
    # first let's fill in the x-edges that are not in tree
    for m = 1:Nx_reducedPhi
        global_ind         = m + ne_x;
        reduced_ind_map[m] = global_ind;
        z[global_ind, m]   = 1.0;
    end
    # then fill in the y-edges that are not in tree
    for m = 1:Ny_reducedPhi
        global_ind = Ne_x + Ney_xyplane + m;
        reduced_ind_map[m+Nx_reducedPhi] = global_ind;
        z[global_ind, m+Nx_reducedPhi]   = 1.0;
    end
    
    # now, the first row of x-edges (these edges are tree branches)
    #for i = 2:ne_x-1
    for i = 1:ne_x
        # first consider the edges in the first xy plane
        for j = 1:ne_y
            tj = j*ne_x + i - ne_x; # subtract ne_x because the indexing of the non-tree edges is shifted down by ne_x compared to their global index
            z[i,tj] = -1.0;
        end
        # now consider edges in other planes
        for j=1:ne_z
            for k = 1:Ny
                tj = j*Nex_xyplane + (k-1)*ne_x +i - ne_x; # subtract ne_x bc the indexing of the non-tree edges is shifted down by ne_x    
                z[i,tj] = -1.0;
            end
        end
    end
    
    # now take care of the y edges in the first column (in y direction) in the lowest xy plane
    for i=1:ne_y
        global_ind = Ne_x + (i-1)*Nx + 1;
        for j= 1:Nz
            for k = (i+1):Ny
                ind = (j-1)*Nex_xyplane + (k-1)*ne_x + 1 - ne_x;
                z[global_ind, ind] = ly/lx;
            end
            if (j > 1)
                indy_global  = Ne_x + (j-1)*Ney_xyplane + (i-1)*Nx + 1;
                indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
                z[global_ind, indy_reduced] = -1.0;
            end
        end
    end
    
    # y edges in the last column (in y direction) in the lowest xy plane 
    for i=1:ne_y
        global_ind = Ne_x + i*Nx; 
        for j=1:Nz
            for k = i+1:Ny
                ind = (j-1)*Nex_xyplane + k*ne_x - ne_x;
                z[global_ind, ind] = -ly/lx;
            end
            if j > 1
                indy_global  = Ne_x + (j-1)*Ney_xyplane + i*Nx;
                indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
                z[global_ind, indy_reduced] = -1.0;
            end
        end
    end
    
    # other y edges in the lowest xy plane
    for i = 2:Nx-1
        for j = 1:ne_y
            global_ind = Ne_x + (j-1)*Nx + i;
            for m = 1:Nz
                for k = j+1:Ny
                    ind1 = (m-1)*Nex_xyplane + (k-1)*ne_x + (i-1) - ne_x;
                    ind2 = ind1 + 1;
                    z[global_ind, ind1] = -ly/lx;
                    z[global_ind, ind2] =  ly/lx;
                end
                if m > 1
                    indy_global  = Ne_x + (m-1)*Ney_xyplane + (j-1)*Nx + i;
                    indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
                    z[global_ind, indy_reduced] = -1.0;
                end
            end
        end
    end
    
    ### Now consider the z edges in the leftmost yz plane (x = xmin plane)
    # first consider the y=ymin column of z edges
    for i = 1:ne_z
        global_ind = Ne_x + Ne_y + (i-1)*Nv_xyplane + 1;
        for j = i+1:Nz
            indx_global  = (j-1)*Nex_xyplane + 1;
            indy_global  = Ne_x + (j-1)*Ney_xyplane + 1;
            indx_reduced = findfirst(x->x== indx_global, reduced_ind_map);
            indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
            z[global_ind, indx_reduced] = lz/lx;
            z[global_ind, indy_reduced] = lz/ly;
        end
    end
    # consider the y=ymax column of z edges
    for i = 1:ne_z
        global_ind = Ne_x + Ne_y + (i-1)*Nv_xyplane + (Ny-1)*Nx + 1;
        for j = i+1:Nz
            indx_global = (j-1)*Nex_xyplane + (Ny-1)*ne_x + 1;
            indy_global = Ne_x + (j-1)*Ney_xyplane + (ne_y-1)*Nx + 1;
            indx_reduced = findfirst(x->x== indx_global, reduced_ind_map);
            indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
            z[global_ind, indx_reduced] =  lz/lx;
            z[global_ind, indy_reduced] = -lz/ly;
        end
    end
    # other z edges in this yz plane
    for i = 2:Ny-1
        for j = 1:ne_z
            global_ind = Ne_x + Ne_y + (j-1)*Nv_xyplane + (i-1)*Nx + 1;
            for k = j+1:Nz
                indx_global   = (k-1)*Nex_xyplane + (i-1)*ne_x + 1;
                indy1_global  = Ne_x + (k-1)*Ney_xyplane + (i-2)*Nx + 1;
                indy2_global  = indy1_global + Nx;
                indx_reduced  = findfirst(x->x== indx_global, reduced_ind_map);
                indy1_reduced = findfirst(x->x== indy1_global, reduced_ind_map);
                indy2_reduced = findfirst(x->x== indy2_global, reduced_ind_map);
                z[global_ind, indx_reduced]  =  lz/lx;
                z[global_ind, indy1_reduced] = -lz/ly;
                z[global_ind, indy2_reduced] =  lz/ly;
            end
        end
    end
    
    ### Now consider the z edges in the rightmost yz plane (x = xmax plane)
    # consider the y=ymin column of z edges
    for i = 1:ne_z
        global_ind = Ne_x + Ne_y + (i-1)*Nv_xyplane + Nx;
        for j = i+1:Nz
            indx_global = (j-1)*Nex_xyplane + ne_x;
            indy_global = Ne_x + (j-1)*Ney_xyplane + Nx;
            indx_reduced = findfirst(x->x== indx_global, reduced_ind_map);
            indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
            z[global_ind, indx_reduced] = -lz/lx;
            z[global_ind, indy_reduced] =  lz/ly;
        end
    end
    # consider the y=ymax column of z edges
    for i = 1:ne_z
        global_ind = Ne_x + Ne_y + i*Nv_xyplane;
        for j = i+1:Nz
            indx_global = j*Nex_xyplane;
            indy_global = Ne_x + j*Ney_xyplane;
            indx_reduced = findfirst(x->x== indx_global, reduced_ind_map);
            indy_reduced = findfirst(x->x== indy_global, reduced_ind_map);
            z[global_ind, indx_reduced] = -lz/lx;
            z[global_ind, indy_reduced] = -lz/ly;
        end
    end
    # consider the other z edges in this yz plane
    for i = 2:Ny-1
        for j = 1:ne_z
            global_ind = Ne_x + Ne_y + (j-1)*Nv_xyplane + i*Nx;
            for k = j+1:Nz
                indx_global   = (k-1)*Nex_xyplane + i*ne_x;
                indy1_global  = Ne_x + (k-1)*Ney_xyplane + (i-1)*Nx;
                indy2_global  = indy1_global + Nx;
                indx_reduced  = findfirst(x->x== indx_global, reduced_ind_map);
                indy1_reduced = findfirst(x->x== indy1_global, reduced_ind_map);
                indy2_reduced = findfirst(x->x== indy2_global, reduced_ind_map);
                z[global_ind, indx_reduced]  = -lz/lx;
                z[global_ind, indy1_reduced] = -lz/ly;
                z[global_ind, indy2_reduced] =  lz/ly;
            end
        end
    end
    
    ### Now consider z edges in other yz planes
    # consider the lowest xz plane of zedges (y = ymin plane)
    for i = 2:Nx-1
        for j = 1:ne_z
            global_ind = Ne_x + Ne_y + (j-1)*Nv_xyplane + i;
            for k = j+1:Nz
                indx1_global = (k-1)*Nex_xyplane + i-1;
                indx2_global = indx1_global + 1;
                indy_global  = Ne_x + (k-1)*Ney_xyplane + i;
                indx1_reduced = findfirst(x->x== indx1_global, reduced_ind_map);
                indx2_reduced = findfirst(x->x== indx2_global, reduced_ind_map);
                indy_reduced  = findfirst(x->x== indy_global, reduced_ind_map);
                z[global_ind, indx1_reduced] = -lz/lx;
                z[global_ind, indx2_reduced] =  lz/lx;
                z[global_ind, indy_reduced]  =  lz/ly;
            end
        end
    end
    # consider the top xz plane of zedges (y = ymax plane)
    for i = 2:Nx-1
        for j = 1:ne_z
            global_ind = Ne_x + Ne_y + (j-1)*Nv_xyplane + Nx*ne_y + i;
            for k = j+1:Nz
                indx1_global = (k-1)*Nex_xyplane + ne_y*ne_x + i-1;
                indx2_global = indx1_global + 1;
                indy_global  = Ne_x + (k-1)*Ney_xyplane + (ne_y-1)*Nx + i;
                indx1_reduced = findfirst(x->x== indx1_global, reduced_ind_map);
                indx2_reduced = findfirst(x->x== indx2_global, reduced_ind_map);
                indy_reduced  = findfirst(x->x== indy_global, reduced_ind_map);
                z[global_ind, indx1_reduced] = -lz/lx;
                z[global_ind, indx2_reduced] =  lz/lx;
                z[global_ind, indy_reduced]  = -lz/ly;
            end
        end
    end
    # consider the remaining z edges
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 1:ne_z
                global_ind = Ne_x + Ne_y + (k-1)*Nv_xyplane + (j-1)*Nx + i;
                for m = k+1:Nz
                    indx1_global = (m-1)*Nex_xyplane + (j-1)*ne_x + i-1;
                    indx2_global = indx1_global + 1;
                    indy1_global = Ne_x + (m-1)*Ney_xyplane + (j-2)*Nx + i;
                    indy2_global = indy1_global + Nx;
                    indx1_reduced = findfirst(x->x== indx1_global, reduced_ind_map);
                    indx2_reduced = findfirst(x->x== indx2_global, reduced_ind_map);
                    indy1_reduced = findfirst(x->x== indy1_global, reduced_ind_map);
                    indy2_reduced = findfirst(x->x== indy2_global, reduced_ind_map);
                    z[global_ind, indx1_reduced] = -lz/lx;
                    z[global_ind, indx2_reduced] =  lz/lx;
                    z[global_ind, indy1_reduced] = -lz/ly;
                    z[global_ind, indy2_reduced] =  lz/ly;
                end
            end
        end
    end
    
    return z;
end
