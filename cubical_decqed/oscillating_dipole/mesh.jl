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

function regionsort3D_cav(e::Array{Int,2},v::Array{Float64,2},Nx::Int64,Ny::Int64,Nv::Int64,ne_x::Int64,ne_y::Int64,Ne_x::Int64,Ne_y::Int64, lx::Float64,ly::Float64, ne_airx::Int64, ne_airy::Int64, ne_airz::Int64, Ne_airx::Int64, Ne_airy::Int64,  Ne_airz::Int64, Ne_air::Int64, Ne_airx_xyplane::Int64, Ne_airy_xyplane::Int64, Ne_airz_xyplane::Int64, Ne_airx_xzplane::Int64, Ne_airz_xzplane::Int64, Ne_airy_yzplane::Int64, Ne_airz_yzplane::Int64, Ne_scx::Int64, Ne_scy::Int64, Ne_scz::Int64, Ne_sc::Int64, Nv_xyplane::Int64, Nex_xyplane::Int64, Ney_xyplane::Int64, ne_scx_seg::Int64, ne_scy_seg::Int64,ne_scz_seg::Int64)

    e_air = Array{Int,1}(undef, Ne_air);
    e_sc  = Array{Int,1}(undef, Ne_sc);

    # find the edges in air
    for i=1:Ne_airx  # find all the x-edges in air, including the edges at interface with sc
        ie_air_xyplane = (i-1)%Ne_airx_xyplane + 1;
        airplane_ind   = div(i-1,Ne_airx_xyplane) + 1;
        
        row_num       = div(ie_air_xyplane-1,ne_airx) + 1; # row in xy superconducting plane
        col_num       = (ie_air_xyplane-1)%ne_airx+1;  # column in xy superconducting plane
        e_air[i] = (ne_scz_seg + airplane_ind -1)*Nex_xyplane + (row_num + ne_scy_seg -1)*ne_x + (col_num + ne_scx_seg);
    end

    for i=1:Ne_airy  # find all the y-edges in air, including the edges at interface with sc
        ie_air_xyplane = (i-1)%Ne_airy_xyplane + 1;
        airplane_ind   = div(i-1,Ne_airy_xyplane) + 1;
        
        row_num       = div(ie_air_xyplane-1,ne_airx+1) + 1;
        col_num       = (ie_air_xyplane-1)%(ne_airx+1) + 1;
        e_air[i+Ne_airx] = Ne_x + (ne_scz_seg + airplane_ind -1)*Ney_xyplane + (row_num + ne_scy_seg -1)*Nx + (col_num + ne_scx_seg);
    end

    for i=1:Ne_airz
        ie_air_xyplane = (i-1)%Ne_airz_xyplane + 1;
        airplane_ind   = div(i-1,Ne_airz_xyplane) + 1;
        
        row_num       = div(ie_air_xyplane-1,ne_airx+1) + 1;
        col_num       = (ie_air_xyplane-1)%(ne_airx+1) + 1;
        e_air[i+Ne_airx+Ne_airy] = Ne_x + Ne_y + (ne_scz_seg + airplane_ind -1)*Nv_xyplane + (row_num + ne_scy_seg -1)*Nx + (col_num + ne_scx_seg);
    end
    
    # find the edges in superconductor
    e_sc_count = 1;
    for i=1:Ne_x  # find all the x-edges in superconductor
        if !(i in e_air)
            e_sc[e_sc_count] = i;
            e_sc_count += 1;
        end
    end

    for i=1:Ne_y  # find all the y-edges in superconductor
        if !(i+Ne_x in e_air)
            e_sc[e_sc_count] = i + Ne_x;
            e_sc_count += 1 ;
        end
    end
    
    for i=1:Ne_z  # find all the z-edges in superconductor
        if !(i+Ne_x+Ne_y in e_air)
            e_sc[e_sc_count] = i + Ne_x + Ne_y;
            e_sc_count += 1 ;
        end
    end
    # Also include into e_sc the edges at the interface with air
    # bottom and top xy interfaces
    for i=1:Ne_airx_xyplane
        row_num = div(i-1,ne_airx) + 1; 
        col_num = (i-1)%ne_airx+1; 
        e_sc[e_sc_count]   = ne_scz_seg*Nex_xyplane + (row_num + ne_scy_seg -1)*ne_x + (col_num + ne_scx_seg);
        e_sc[e_sc_count+1] = (ne_scz_seg+ne_airz)*Nex_xyplane + (row_num + ne_scy_seg -1)*ne_x + (col_num + ne_scx_seg);
        e_sc_count += 2;
    end
    for i=1:Ne_airy_xyplane
        row_num = div(i-1,ne_airx+1) + 1;
        col_num = (i-1)%(ne_airx+1) + 1;
        e_sc[e_sc_count]   = Ne_x + ne_scz_seg*Ney_xyplane + (row_num + ne_scy_seg -1)*Nx + (col_num + ne_scx_seg);
        e_sc[e_sc_count+1] = Ne_x + (ne_scz_seg+ne_airz)*Ney_xyplane + (row_num + ne_scy_seg -1)*Nx + (col_num + ne_scx_seg);
        e_sc_count += 2;
    end
    #  xz interfaces
    for i=1:(Ne_airx_xzplane-2*ne_airx) # only count from 2nd to 2nd-to-last row
        airplane_ind = div(i-1,ne_airx) + 2; # starts from the second row, since the 1st row is already counted in xy plane
        col_num      = (i-1)%ne_airx + 1;
        e_sc[e_sc_count]   = (ne_scz_seg + airplane_ind -1)*Nex_xyplane + ne_scy_seg*ne_x + (col_num + ne_scx_seg);
        e_sc[e_sc_count+1] = (ne_scz_seg + airplane_ind -1)*Nex_xyplane + (ne_airy + ne_scy_seg)*ne_x + (col_num + ne_scx_seg);
        e_sc_count += 2;
    end
    for i=1:Ne_airz_xzplane
        airplane_ind = div(i-1,ne_airx+1) + 1;
        col_num      = (i-1)%(ne_airx+1) + 1;
        e_sc[e_sc_count]   = Ne_x + Ne_y + (ne_scz_seg + airplane_ind -1)*Nv_xyplane + ne_scy_seg*Nx + (col_num + ne_scx_seg);
        e_sc[e_sc_count+1] = Ne_x + Ne_y + (ne_scz_seg + airplane_ind -1)*Nv_xyplane + (ne_airy + ne_scy_seg)*Nx + (col_num + ne_scx_seg);
        e_sc_count += 2;
    end
    #  yz interfaces
    for i=1:(Ne_airy_yzplane-2*ne_airy) # only count from 2nd to 2nd-to-last row
        airplane_ind = div(i-1,ne_airy) + 2; # starts from the second row, since the 1st row is already counted in xy plane
        row_num      = (i-1)%ne_airy + 1;
        e_sc[e_sc_count]   = Ne_x + (ne_scz_seg + airplane_ind -1)*Ney_xyplane + (row_num + ne_scy_seg -1)*Nx + (1 + ne_scx_seg);
        e_sc[e_sc_count+1] = Ne_x + (ne_scz_seg + airplane_ind -1)*Ney_xyplane + (row_num + ne_scy_seg -1)*Nx + (1 + ne_airx + ne_scx_seg);
        e_sc_count += 2;
    end
    for i=1:(Ne_airz_yzplane-2*ne_airz) # only count from 2nd to 2nd-to-last column
        airplane_ind = div(i-1,ne_airy-1) + 1; # starts from the second column, since the 1st column is already counted in xz plane
        row_num      = (i-1)%(ne_airy-1) + 2;  # starts from the second column, since the 1st column is already counted in xz plane
        e_sc[e_sc_count]   = Ne_x + Ne_y + (ne_scz_seg + airplane_ind -1)*Nv_xyplane + (row_num + ne_scy_seg -1)*Nx + (1 + ne_scx_seg);
        e_sc[e_sc_count+1] = Ne_x + Ne_y + (ne_scz_seg + airplane_ind -1)*Nv_xyplane + (row_num + ne_scy_seg -1)*Nx + (1 + ne_airx + ne_scx_seg);
        e_sc_count += 2;
    end
    
    
    return e_sc, e_air;
end


function materials_3Dcav(Ne::Int64, e_sc::Array{Int,1}, lambda::Float64)
    invLambda2 = Array{Float64,1}(undef, Ne);
    for i = 1:Ne
        if (i in e_sc)
            invLambda2[i] = 1/lambda^2;
        else
            invLambda2[i] = 0.0;
        end
    end
    return invLambda2;
end
