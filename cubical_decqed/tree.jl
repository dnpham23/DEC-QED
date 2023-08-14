function stree3D(ne_x::Int64,Nbranch::Int64,Ney_xyplane::Int64,Ne_x::Int64,Ne_y::Int64, Ne_z::Int64)
    #tree = Array{Int,2}(undef, Nbranch, 1);
    tree = Array{Int,1}(undef, Nbranch);
    # first, take the x edges lying on the y=ymin, z=zmin line of the domain boundary
    for i = 1:ne_x 
        tree[i] = i;
    end
    
    # take the y edges lying on the lowest xy plane (the z=zmin plane) of the domain boundary
    for i = 1:Ney_xyplane
        tree[i+ne_x] = Ne_x + i;
    end
    
    # take all the z edges in the domain
    for i = 1:Ne_z
        tree[i+ne_x+Ney_xyplane] = Ne_x + Ne_y +i;
    end
    
    return tree;
end