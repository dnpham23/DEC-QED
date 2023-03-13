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
function eboundary3D(vbound::Array{Int,1},e::Array{Int,2},ne_x::Int64,ne_y::Int64,ne_z::Int64,Ne::Int64,Nx::Int64, Ny::Int64, Nz::Int64, Ne_x::Int64, Ne_y::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64)
    Nebound = 2*Nex_xyplane + 2*ne_x*(Nz-2) + 2*Ney_xyplane + 2*ne_y*(Nz-2) + 2*ne_z*Nx + 2*ne_z*(Ny-2);
    Nebound_all = Nebound + 2*(ne_x-1)*(ne_y-1) + 2*(ne_x-1)*(ne_z-1) + 2*(ne_y-1)*(ne_z-1);
    
    ebound     = Array{Int,1}(undef, Nebound);    
    bcount = 1;
    for i = 1:Ne
        if (e[i,1] in vbound)&&(e[i,2] in vbound)
            ebound[bcount] = i;
            bcount += 1;
        end
    end
    
    ebound_all = Array{Int,1}(undef, Nebound_all);
    ebound_all[1:length(ebound)] = ebound;
    ebound_all_cnt = length(ebound) + 1;

    for xind=1:Ne_x
        ie_xyplane = (xind-1)%Nex_xyplane + 1;
        plane_ind  = div(xind-1,Nex_xyplane) + 1;
        rownum     = div(ie_xyplane-1,ne_x) + 1; # column in the xy plane
        colnum     = (ie_xyplane-1)%ne_x + 1;    # row in the xy plane
        if ((colnum==1)||(colnum==ne_x))&&(rownum!=1)&&(rownum!=Ny)&&(plane_ind!=1)&&(plane_ind!=Nz)
            ebound_all[ebound_all_cnt] = xind;
            ebound_all_cnt += 1;
        end
    end
    
    for yind=1:Ne_y
        real_yind = yind + Ne_x;
        ie_xyplane = (yind-1)%Ney_xyplane + 1;
        plane_ind  = div(yind-1,Ney_xyplane) + 1;
        rownum     = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
        colnum     = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
        if ((rownum==1)||(rownum==ne_y))&&(colnum!=1)&&(colnum!=Nx)&&(plane_ind!=1)&&(plane_ind!=Nz)
            ebound_all[ebound_all_cnt] = real_yind;
            ebound_all_cnt += 1; 
        end
    end
    
    for zind=1:Ne_z
        real_zind = zind + Ne_x + Ne_y;
        ie_xyplane = (zind-1)%Nv_xyplane + 1;
        plane_ind  = div(zind-1,Nv_xyplane) + 1;
        rownum     = div(ie_xyplane-1,Nx) + 1; # column in the xy plane
        colnum     = (ie_xyplane-1)%Nx + 1;    # row in the xy plane
        if ((plane_ind==1)||(plane_ind==ne_z))&&(colnum!=1)&&(colnum!=Nx)&&(rownum!=1)&&(rownum!=Ny)
            ebound_all[ebound_all_cnt] = real_zind;
            ebound_all_cnt += 1;             
        end
    end
    
    return ebound, ebound_all;
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

function regionsort3D_scring(e::Array{Int,2},v::Array{Float64,2},Nx::Int64,Ny::Int64,Nv::Int64,ne_x::Int64,ne_y::Int64,Ne_x::Int64,Ne_y::Int64, lx::Float64,ly::Float64, Ne_scx::Int64, Ne_scy::Int64,Ne_scz::Int64, Ne_sc::Int64, Nv_xyplane::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64, ne_airx_seg::Int64, ne_airy_seg::Int64,ne_airz_seg::Int64, ne_scx_horblock::Int64, ne_scx_verblock::Int64, ne_scy_verblock::Int64, ne_scy_horblock::Int64, nv_sc_horblock::Int64, nv_sc_verblock::Int64, Ne_scx_xyplane::Int64, Ne_scy_xyplane::Int64, Nv_sc_xyplane::Int64)

    #e_air = Array{Int,1}(undef, Ne_air);
    e_sc  = Array{Int,1}(undef, Ne_sc);

    # find the edges in sc
    for iz = 1:(ne_scz+1)
        sc_xbase = (iz-1)*Ne_scx_xyplane;
        for i=1:ne_scx_horblock
            rownum = div(i-1,ne_scx_len)+1;
            colnum = (i-1)%ne_scx_len+1;
            e_sc[sc_xbase+i] = (iz+ne_airz_seg-1)*Nex_xyplane + (rownum-1+ne_airy_seg)*ne_x + colnum + ne_airx_seg;
            e_sc[sc_xbase+ne_scx_horblock+i] = (iz+ne_airz_seg-1)*Nex_xyplane + (rownum-1+ne_airy_seg+ne_scy_seg+ne_scy_inlen)*ne_x + colnum + ne_airx_seg;
        end
        
        for i=1:ne_scx_verblock
            rownum = div(i-1,ne_scx_seg)+1;
            colnum = (i-1)%ne_scx_seg + 1;
            e_sc[sc_xbase+2*ne_scx_horblock+i] = (iz+ne_airz_seg-1)*Nex_xyplane + (rownum-1+ne_scy_seg+1+ne_airy_seg)*ne_x + colnum + ne_airx_seg;
            e_sc[sc_xbase+2*ne_scx_horblock+ne_scx_verblock+i] = (iz+ne_airz_seg-1)*Nex_xyplane + (rownum-1+ne_scy_seg+1+ne_airy_seg)*ne_x + colnum + ne_airx_seg + ne_scx_seg + ne_scx_inlen;
        end
        
        sc_ybase = (iz-1)*Ne_scy_xyplane;
        for i=1:ne_scy_verblock
            rownum = div(i-1,ne_scx_seg+1)+1;
            colnum = (i-1)%(ne_scx_seg+1)+1;
            e_sc[Ne_scx+sc_ybase+i] = Ne_x + (iz+ne_airz_seg-1)*Ney_xyplane + (rownum+ne_airy_seg-1)*Nx + colnum + ne_airx_seg;
            e_sc[Ne_scx+sc_ybase+ne_scy_verblock+i] = Ne_x + (iz+ne_airz_seg-1)*Ney_xyplane + (rownum+ne_airy_seg-1)*Nx + colnum + ne_airx_seg + ne_scx_seg + ne_scx_inlen;
        end
        
        for i=1:ne_scy_horblock
            rownum = div(i-1,ne_scx_inlen-1)+1;
            colnum = (i-1)%(ne_scx_inlen-1)+1;
            e_sc[Ne_scx+sc_ybase+2*ne_scy_verblock+i] = Ne_x + (iz+ne_airz_seg-1)*Ney_xyplane + (rownum-1+ne_airy_seg)*Nx + colnum + ne_airx_seg + ne_scx_seg + 1;
            e_sc[Ne_scx+sc_ybase+2*ne_scy_verblock+ne_scy_horblock+i] = Ne_x + (iz+ne_airz_seg-1)*Ney_xyplane + (rownum-1+ne_airy_seg+ne_scy_seg+ne_scy_inlen)*Nx + colnum + ne_airx_seg + ne_scx_seg + 1;
        end
    end
    
    for iz = 1:ne_scz
        sc_zbase = (iz-1)*Nv_sc_xyplane;
        for i = 1:nv_sc_horblock
            rownum = div(i-1,ne_scx_len+1)+1;
            colnum = (i-1)%(ne_scx_len+1) + 1;
            e_sc[Ne_scx+Ne_scy+sc_zbase+i] = Ne_x + Ne_y + (iz+ne_airz_seg-1)*Nv_xyplane + (rownum-1+ne_airy_seg)*Nx + colnum + ne_airx_seg;
            e_sc[Ne_scx+Ne_scy+sc_zbase+nv_sc_horblock+i] = Ne_x + Ne_y + (iz+ne_airz_seg-1)*Nv_xyplane + (rownum-1+ne_airy_seg+ne_scy_seg+ne_scy_inlen)*Nx + colnum + ne_airx_seg;
        end
        
        for i = 1:nv_sc_verblock
            rownum = div(i-1,ne_scx_seg+1)+1;
            colnum = (i-1)%(ne_scx_seg+1)+1;
            e_sc[Ne_scx+Ne_scy+sc_zbase+2*nv_sc_horblock+i] = Ne_x + Ne_y + (iz+ne_airz_seg-1)*Nv_xyplane + (rownum-1+ne_scy_seg+1+ne_airy_seg)*Nx + colnum + ne_airx_seg;
            e_sc[Ne_scx+Ne_scy+sc_zbase+2*nv_sc_horblock+nv_sc_verblock+i] = Ne_x + Ne_y + (iz+ne_airz_seg-1)*Nv_xyplane + (rownum-1+ne_scy_seg+1+ne_airy_seg)*Nx + colnum + ne_airx_seg + ne_scx_seg + ne_scx_inlen;
        end
    end
    
    return e_sc;
end


function materials_scringJJ(Ne::Int64, e_sc::Array{Int,1}, e_jj::Array{Int,1}, lambda::Float64)
    invLambda2 = Array{Float64,1}(undef, Ne);
    for i = 1:Ne
        if (i in e_sc)&&(!(i in e_jj))
            invLambda2[i] = 1/lambda^2;
        else
            invLambda2[i] = 0.0;
        end
    end
    return invLambda2;
end
