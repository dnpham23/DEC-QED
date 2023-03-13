function phip_grid_distribute(phip_next_x::Array{Float64,1},phip_next_y::Array{Float64,1},ne_scx_seg::Int64,ne_scy_seg::Int64,Nx::Int64, Ne_x::Int64, Ne_y::Int64,e_scx::Array{Int,1},e_scy::Array{Int,1}, Ne_scx::Int64, Ne_scy::Int64,ne_x::Int64, ne_y::Int64, Ny::Int64)
    phipx_vec_next = zeros(Ne_x,1);
    phipy_vec_next = zeros(Ne_y,1);
    
    # fill in the phipx grid
    for i=1:Ne_scx
        real_ind = e_scx[i];
        phipx_vec_next[real_ind] = phip_next_x[i];
    end
    
    for i=1:Ne_scy
        real_ind = e_scy[i];
        phipy_vec_next[real_ind-Ne_x] = phip_next_y[i];
    end
    phipx_grid_next = reshape(phipx_vec_next,ne_x,Ny);
    phipy_grid_next = reshape(phipy_vec_next,Nx,ne_y);
    
    return phipx_grid_next, phipy_grid_next;
end