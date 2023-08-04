# Generate a rectangular cuboid mesh
function Mesh3Dcube(xmax::Float64,xmin::Float64,ymax::Float64,ymin::Float64,zmax::Float64,zmin::Float64,Nx::Int64,Ny::Int64,Nz::Int64, ne_x::Int64,ne_y::Int64,ne_z::Int64,Ne_x::Int64,Ne_y::Int64, Ne_z::Int64, Ne::Int64,Nv::Int64,lx::Float64,ly::Float64,lz::Float64, Nv_xyplane::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64)

    v = Array{Float64,2}(undef, Nv, 3);
    e = Array{Int,2}(undef, Ne, 2);
    for i = 1:Nv
        iv_xyplane = (i-1)%Nv_xyplane + 1;
        plane_ind  = div(i-1,Nv_xyplane) + 1;
        v[i,1] = xmin + ((iv_xyplane-1)%Nx)*lx;
        v[i,2] = ymin + div(iv_xyplane-1,Nx)*ly;
        v[i,3] = zmin + (plane_ind-1)*lz;
    end

    for i = 1:Ne_x
        ie_xyplane = (i-1)%Nex_xyplane + 1;
        plane_ind  = div(i-1,Nex_xyplane) + 1;
        etemp = div(ie_xyplane-1,ne_x)*Nx + (ie_xyplane-1)%ne_x + 1 + (plane_ind-1)*Nv_xyplane;
        e[i,1] = etemp;
        e[i,2] = etemp+1;
    end

    for i = 1:Ne_y
        ie_xyplane = (i-1)%Ney_xyplane + 1;
        plane_ind  = div(i-1,Ney_xyplane) + 1;
        e[i+Ne_x,1] = ie_xyplane + (plane_ind-1)*Nv_xyplane;
        e[i+Ne_x,2] = ie_xyplane + (plane_ind-1)*Nv_xyplane + Nx;
    end
    
    for i = 1:Ne_z
        ie_xyplane = (i-1)%Nv_xyplane + 1;
        plane_ind  = div(i-1,Nv_xyplane) + 1;
        e[i+Ne_x+Ne_y,1] = ie_xyplane + (plane_ind-1)*Nv_xyplane;
        e[i+Ne_x+Ne_y,2] = ie_xyplane + plane_ind*Nv_xyplane;
    end
    
    return e,v;
end

# find the boundary vertices
function vboundary3D(v::Array{Float64,2},xmax::Float64,xmin::Float64,ymax::Float64,ymin::Float64,zmax::Float64,zmin::Float64, Nv::Int64,Nx::Int64,Ny::Int64,Nv_xyplane::Int64)
    xtol = (xmax-xmin)*1e-5;
    ytol = (ymax-ymin)*1e-5;
    ztol = (zmax-zmin)*1e-5;
    Nvbound = 2*Nv_xyplane + 2*Nx*(Nz-2) + 2*(Ny-2)*(Nz-2);
    vbound = Array{Int,1}(undef, Nvbound);
    bcount = 1;
    for i = 1:Nv
        if (abs(v[i,1]-xmin)<xtol)||(abs(v[i,1]-xmax)<xtol)||(abs(v[i,2]-ymin)<ytol)||(abs(v[i,2]-ymax)<ytol)||(abs(v[i,3]-zmin)<ztol)||(abs(v[i,3]-zmax)<ztol)
            vbound[bcount] = i;
            bcount +=1;
        end
    end
    return vbound, xtol, ytol, ztol;
end

# find the boundary edges
function eboundary(vbound::Array{Int,1},e::Array{Int,2},ne_x::Int64,ne_y::Int64,ne_z::Int64,Ne::Int64,Ne_sc::Int64, Ne_scx::Int64,Ne_scy::Int64, e_scx::Array{Int,1},e_scy::Array{Int,1},Nx::Int64, Ny::Int64,Ne_x::Int64)
    Nebound = 2*(ne_x+ne_y);
    ebound = Array{Int,1}(undef, Nebound);
    bcount = 1;
    for i = 1:Ne
        if (e[i,1] in vbound)&&(e[i,2] in vbound)
            ebound[bcount] = i;
            bcount += 1;
        end
    end
    
    ebound_all = Array{Int,1}(undef, length(ebound)+ 2*(Nx-2 + Ny-2));
    ebound_all[1:length(ebound)] = ebound;
    ebound_all_cnt = length(ebound) + 1;
    for i = 1:Ne_sc
        if i <= Ne_scx
            real_ind   = e_scx[i];
            rownum     = div(real_ind-1,ne_x) + 1; # global row
            colnum     = (real_ind-1)%ne_x + 1;    # global column
            if ((colnum==1)||(colnum==ne_x))&&(rownum!=1)&&(rownum!=Ny)
                ebound_all[ebound_all_cnt] = real_ind;
                ebound_all_cnt += 1;
            end
        else
            real_ind = e_scy[i-Ne_scx];
            rownum    = div(real_ind-Ne_x-1,Nx) + 1;
            colnum    = (real_ind-Ne_x-1)%Nx + 1;
            if ((rownum==1)||(rownum==ne_y))&&(colnum!=1)&&(colnum!=Nx)
                ebound_all[ebound_all_cnt] = real_ind;
                ebound_all_cnt += 1;
            end
        end
    end
    return ebound, ebound_all;
end

# find the boundary edges
function eboundary3D(vbound::Array{Int,1},e::Array{Int,2},ne_x::Int64,ne_y::Int64,ne_z::Int64,Ne::Int64,Nx::Int64, Ny::Int64, Nz::Int64, Ne_x::Int64, Ne_y::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64)
    Nebound = 2*Nex_xyplane + 2*ne_x*(Nz-2) + 2*Ney_xyplane + 2*ne_y*(Nz-2) + 2*ne_z*Nx + 2*ne_z*(Ny-2);
    ebound = Array{Int,1}(undef, Nebound);
    bcount = 1;
    for i = 1:Ne
        if (e[i,1] in vbound)&&(e[i,2] in vbound)
            ebound[bcount] = i;
            bcount += 1;
        end
    end
    
    return ebound;
end

# Map each vertex to the edges that are attached to it
function vemap3D(e::Array{Int,2},v::Array{Float64,2},vbound::Array{Int,1},Nx::Int64,Ny::Int64,Nz::Int64,Nv::Int64,ne_x::Int64,Ne_x::Int64, Ne_y::Int64, xtol::Float64, ytol::Float64,ztol::Float64,xmin::Float64,xmax::Float64,ymin::Float64,ymax::Float64,zmin::Float64,zmax::Float64, Nex_xyplane::Int64, Ney_xyplane::Int64, Nv_xyplane::Int64)
    v2exmap = Array{Int,2}(undef, Nv, 2);
    v2eymap = Array{Int,2}(undef, Nv, 2);
    v2ezmap = Array{Int,2}(undef, Nv, 2);    
    for i = 1:Nv
        iv_xyplane = (i-1)%Nv_xyplane + 1;
        plane_ind  = div(i-1,Nv_xyplane) + 1;
        vrow       = div(iv_xyplane-1,Nx) + 1; # row in xy plane
        vcol       = (iv_xyplane-1)%Nx + 1;    # column in xy plane
            
        v2exmap[i,1] = ne_x*(vrow-1) + vcol-1 + (plane_ind-1)*Nex_xyplane;    # -x edge
        v2exmap[i,2] = ne_x*(vrow-1) + vcol + (plane_ind-1)*Nex_xyplane;      # +x edge
        v2eymap[i,1] = Ne_x + Nx*(vrow-2) + vcol + (plane_ind-1)*Ney_xyplane; # -y edge
        v2eymap[i,2] = Ne_x + Nx*(vrow-1) + vcol + (plane_ind-1)*Ney_xyplane; # +y edge
        v2ezmap[i,1] = Ne_x + Ne_y + iv_xyplane + (plane_ind-2)*Nv_xyplane;  # -z edge
        v2ezmap[i,2] = Ne_x + Ne_y + iv_xyplane + (plane_ind-1)*Nv_xyplane;  # -z edge
        
        if (abs(v[i,1]-xmin)<xtol)
            v2exmap[i,1] = 0;   # -x edge
        end
        if (abs(v[i,1]-xmax)<xtol)
            v2exmap[i,2] = 0;   # +x edge
        end
        if (abs(v[i,2]-ymin)<ytol)
            v2eymap[i,1] = 0;   # -y edge
        end
        if (abs(v[i,2]-ymax)<ytol)
            v2eymap[i,2] = 0;   # +y edge
        end
        if (abs(v[i,3]-zmin)<ztol)
            v2ezmap[i,1] = 0;   # -z edge
        end
        if (abs(v[i,3]-zmax)<ztol)
            v2ezmap[i,2] = 0;   # +z edge
        end
    end
    return v2exmap, v2eymap, v2ezmap;
end

function regionsort3D(e::Array{Int,2},v::Array{Float64,2},Nx::Int64,Ny::Int64,Nv::Int64,ne_x::Int64,ne_y::Int64,Ne_x::Int64,Ne_y::Int64, lx::Float64,ly::Float64,ne_scx::Int64,ne_scy::Int64,ne_scz::Int64,Ne_scx::Int64,Ne_scy::Int64,Ne_scz::Int64, Ne_sc::Int64, Ne_scx_xyplane::Int64, Ne_scy_xyplane::Int64, Ne_scz_xyplane::Int64, Nv_xyplane::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64, ne_airx_seg::Int64, ne_airy_seg::Int64,ne_airz_seg::Int64)

    #e_air = Array{Int,1}(undef, Ne_air);
    e_sc  = Array{Int,1}(undef, Ne_sc);

    # find the edges in sc
    for i=1:Ne_scx  # find all the x-edges in sc, including the edges at interface with air
        ie_sc_xyplane = (i-1)%Ne_scx_xyplane + 1;
        scplane_ind   = div(i-1,Ne_scx_xyplane) + 1;
        
        row_num       = div(ie_sc_xyplane-1,ne_scx) + 1; # row in xy superconducting plane
        col_num       = (ie_sc_xyplane-1)%ne_scx+1;  # column in xy superconducting plane
        e_sc[i] = (ne_airz_seg + scplane_ind -1)*Nex_xyplane + (row_num + ne_airy_seg -1)*ne_x + (col_num + ne_airx_seg);
    end

    for i=1:Ne_scy  # find all the y-edges in sc, including the edges at interface with air
        ie_sc_xyplane = (i-1)%Ne_scy_xyplane + 1;
        scplane_ind   = div(i-1,Ne_scy_xyplane) + 1;
        
        row_num       = div(ie_sc_xyplane-1,ne_scx+1) + 1;
        col_num       = (ie_sc_xyplane-1)%(ne_scx+1) + 1;
        e_sc[i+Ne_scx] = Ne_x + (ne_airz_seg + scplane_ind -1)*Ney_xyplane + (row_num + ne_airy_seg -1)*Nx + (col_num + ne_airx_seg);
    end

    for i=1:Ne_scz
        ie_sc_xyplane = (i-1)%Ne_scz_xyplane + 1;
        scplane_ind   = div(i-1,Ne_scz_xyplane) + 1;
        
        row_num       = div(ie_sc_xyplane-1,ne_scx+1) + 1;
        col_num       = (ie_sc_xyplane-1)%(ne_scx+1) + 1;
        e_sc[i+Ne_scx+Ne_scy] = Ne_x + Ne_y + (ne_airz_seg + scplane_ind -1)*Nv_xyplane + (row_num + ne_airy_seg -1)*Nx + (col_num + ne_airx_seg);
    end
    
    return e_sc;
end

function materials(Ne::Int64, e_sc::Array{Int,1}, mu_sc::Float64,eps_sc::Float64,mu_air::Float64,eps_air::Float64,lambda::Float64)
    #eps_list   = Array{Float64,1}(undef, Ne);
    #mu_list    = Array{Float64,1}(undef, Ne);
    invLambda2 = Array{Float64,1}(undef, Ne);
    for i = 1:Ne
        if (i in e_sc)
            #eps_list[i]   = eps_sc;
            #mu_list[i]    = mu_sc;
            invLambda2[i] = 1/lambda^2;
        else
            #eps_list[i]   = eps_air;
            #mu_list[i]    = mu_air;
            invLambda2[i] = 0.0;
        end
    end
    return invLambda2;
end
