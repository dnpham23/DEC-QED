"""
    vboundary_sphere(r1::Float64,r2::Float64,r3::Float64,xo::Float64,yo::Float64,zo::Float64, primalmesh::Primalmeshstruct, tol::Float64)

    find boundary nodes
"""
function vboundary_sphere(r1::Float64,r2::Float64,r3::Float64,xo::Float64,yo::Float64,zo::Float64, primalmesh::Primalmeshstruct, tol::Float64)
    
    vbound  = []
    vbound1 = []
    vbound2 = []
    vbound3 = []
    Nvbound_all = 0
    
    for k in keys(primalmesh.nodedict)
        x = primalmesh.nodedict[k].coords[1]
        y = primalmesh.nodedict[k].coords[2]
        z = primalmesh.nodedict[k].coords[3]
        r = sqrt((x-xo)^2 + (y-yo)^2 + (z-zo)^2)
        if (abs(r-r1)/r1<=tol)
            push!(vbound1, k)
            Nvbound_all += 1
        elseif (abs(r-r2)/r2<=tol)
            push!(vbound2, k)
            Nvbound_all += 1
        elseif (abs(r-r3)/r3<=tol)
            push!(vbound3, k)
            Nvbound_all += 1
        end
    end
   
    vbound = [vbound1; vbound2; vbound3]
    return vbound, vbound1, vbound2, vbound3, Nvbound_all
end

"""
    eboundary_sphere(vbound1::Array{Any,1}, vbound2::Array{Any,1}, vbound3::Array{Any,1}, primalmesh::Primalmeshstruct, r1::Float64, r2::Float64, r3::Float64, xo::Float64, yo::Float64, zo::Float64, tol::Float64)

    find boundary edges
"""
function eboundary_sphere(vbound1::Array{Any,1}, vbound2::Array{Any,1}, vbound3::Array{Any,1}, primalmesh::Primalmeshstruct, r1::Float64, r2::Float64, r3::Float64, xo::Float64, yo::Float64, zo::Float64, tol::Float64)
    # define array of edges on the boundary surfaces
    ebound1 = [] 
    ebound2 = []
    ebound3 = []
    # define the array of radial edges that connect two consecutive boundary surfaces
    ebound_r12 = []
    ebound_r23 = [] 
    
    #r12_ratio = (r2-r1)/r2
    #r23_ratio = (r3-r2)/r3
    r12       = r2-r1
    r23       = r3-r2

    for k in keys(primalmesh.edgedict)
        v1 = k[1]
        v2 = k[2]
        x1 = primalmesh.nodedict[v1].coords[1]
        y1 = primalmesh.nodedict[v1].coords[2]
        z1 = primalmesh.nodedict[v1].coords[3]
        x2 = primalmesh.nodedict[v2].coords[1]
        y2 = primalmesh.nodedict[v2].coords[2]
        z2 = primalmesh.nodedict[v2].coords[3]
        dr = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)

        if (v1 in vbound1)&&(v2 in vbound1)
            push!(ebound1, k)
        elseif (v1 in vbound2)&&(v2 in vbound2)
            push!(ebound2, k)
        elseif (v1 in vbound3)&&(v2 in vbound3)
            push!(ebound3, k)
	###elseif ( ((v1 in vbound1)&&(v2 in vbound2))||((v1 in vbound2)&&(v2 in vbound1)) )&&(abs(abs((x1-x2)/x2)-r12_ratio)<=tol)&&(abs(abs((y1-y2)/y2)-r12_ratio)<=tol)&&(abs(abs((z1-z2)/z2)-r12_ratio)<=tol)
	elseif ( ((v1 in vbound1)&&(v2 in vbound2))||((v1 in vbound2)&&(v2 in vbound1)) )&&(abs(dr-r12)/r12<=tol)
            push!(ebound_r12, k)
        ###elseif ( ((v1 in vbound2)&&(v2 in vbound3))||((v1 in vbound3)&&(v2 in vbound2)) )&&(abs(abs((x1-x2)/x2)-r23_ratio)<=tol)&&(abs(abs((y1-y2)/y2)-r23_ratio)<=tol)&&(abs(abs((z1-z2)/z2)-r23_ratio)<=tol)
	elseif ( ((v1 in vbound2)&&(v2 in vbound3))||((v1 in vbound3)&&(v2 in vbound2)) )&&(abs(dr-r23)/r23<=tol)
            push!(ebound_r23, k) 
        end
    end
    return ebound1, ebound2, ebound3, ebound_r12, ebound_r23
end

"""
    fboundary_sphere(vbound1::Array{Any,1}, vbound2::Array{Any,1}, vbound3::Array{Any,1}, primalmesh::Primalmeshstruct)

    find the faces on the boundary layers
"""
function fboundary_sphere(vbound1::Array{Any,1}, vbound2::Array{Any,1}, vbound3::Array{Any,1}, primalmesh::Primalmeshstruct)
    
    fbound1 = []
    fbound2 = []
    fbound3 = []
    
    for k in keys(primalmesh.facedict)
        v1 = k[1]
        v2 = k[2]
        v3 = k[3]
        
        if (v1 in vbound1)&&(v2 in vbound1)&&(v3 in vbound1)
            push!(fbound1, k)
        elseif (v1 in vbound2)&&(v2 in vbound2)&&(v3 in vbound2)
            push!(fbound2, k)
        elseif (v1 in vbound3)&&(v2 in vbound3)&&(v3 in vbound3)
            push!(fbound3, k)
        end
    end
    return fbound1, fbound2, fbound3
end

"""
    bnode_add(vbound3::Array{Any,1}, primalmesh::Primalmeshstruct, r1::Float64, r2::Float64, r3::Float64, xo::Float64, yo::Float64, zo::Float64, tol::Float64)

    add a new layer of boundary nodes outside of boundary layer r3
"""
function bnode_add(vbound3::Array{Any,1}, primalmesh::Primalmeshstruct, r1::Float64, r2::Float64, r3::Float64, xo::Float64, yo::Float64, zo::Float64, tol::Float64)

    # dictionary to store the new boundary nodes
    nodedict_extra = Dict{Int, Nodestruct}()
    vbound4 = []
   
    # the number of nodes on each boundary layer
    Nvbound_layer = length(vbound3)
    # the largest original nodeid value
    nodeid_max = findmax(collect(keys(primalmesh.nodedict)))[1]
    # 
    dr = r3-r2
    
    for nv = 1:Nvbound_layer
        nodeid =  nodeid_max + nv
        # id of the corresponding node on boundary r3
        matching_node_id = vbound3[nv]
        # find the coordinates of the coresponding node on boundary r3
        xm = primalmesh.nodedict[matching_node_id].coords[1]
        ym = primalmesh.nodedict[matching_node_id].coords[2]
        zm = primalmesh.nodedict[matching_node_id].coords[3]
        # define the new node
        nodecoords = [xo + (xm-xo)*(r3+dr)/r3, yo + (ym-yo)*(r3+dr)/r3, zo + (zm-zo)*(r3+dr)/r3]
        node = Nodestruct()
        node.id = nodeid
        node.coords = nodecoords
        node.root_entityid = 0
        nodedict_extra[nodeid] = node
        push!(vbound4, nodeid)
    end 
    return nodedict_extra, vbound4, r3+dr
end

"""
    bedge_add(vbound3::Array{Any,1}, vbound4::Array{Any,1}, primalmesh::Primalmeshstruct, nodedict_extra::Dict{Int64,Nodestruct}, r3::Float64, r4::Float64)

    define new radial edges that connects layer 3 with the new outer boundary nodes
"""
function bedge_add(vbound3::Array{Any,1}, vbound4::Array{Any,1}, primalmesh::Primalmeshstruct, nodedict_extra::Dict{Int64,Nodestruct}, r3::Float64, r4::Float64)
    
    # dictionary to store the new radial boundary edges
    edgedict_extra = Dict{SVector{2, Int}, Edgestruct}()
    
    for ne = 1:length(vbound4)
        node1 = vbound3[ne]
        node2 = vbound4[ne]
        edgeid = [node1, node2]
        edge = Edgestruct()
        edge.id = edgeid
        edge.length = r4-r3
        #edge.length = sqrt( (primalmesh.nodedict[node1].coords[1]-nodedict_extra[node2].coords[1])^2 + (primalmesh.nodedict[node1].coords[2]-nodedict_extra[node2].coords[2])^2 + (primalmesh.nodedict[node1].coords[3]-nodedict_extra[node2].coords[3])^2 ) ### this line is to double check that the lengths are indeed r4-r3
        edgedict_extra[edgeid] = edge
    end
    return edgedict_extra
end

# check whether an edge actually exists
#edgekeylist = collect(keys(primalmesh.edgedict))
#temp = []
#for fkey in fbound1
#    edges_list = primalmesh.facedict[fkey].edges
#    if (!(edges_list[1] in edgekeylist))||(!(edges_list[2] in edgekeylist))||!(edges_list[3] in edgekeylist) 
#        push!(temp, fkey)
#    end
#end


# collect all the edges inside the dielectric sphere
function dielectric_edgecollect(primalmesh::Primalmeshstruct, rd::Float64, tol::Float64, xo::Float64, yo::Float64, zo::Float64, edgekeylist)

    edge_dielectric = []
    edge_dielec_ind = []
    for k in keys(primalmesh.edgedict)
        v1 = k[1]
        v2 = k[2]
        x1 = primalmesh.nodedict[v1].coords[1]
        y1 = primalmesh.nodedict[v1].coords[2]
        z1 = primalmesh.nodedict[v1].coords[3]
	rp1 = sqrt((x1-xo)^2+(y1-yo)^2+(z1-zo)^2)
        x2 = primalmesh.nodedict[v2].coords[1]
        y2 = primalmesh.nodedict[v2].coords[2]
        z2 = primalmesh.nodedict[v2].coords[3]
        rp2 = sqrt((x2-xo)^2+(y2-yo)^2+(z2-zo)^2)

	if ( (rp1<=rd)||(abs(rp1-rd)/rd<=tol) )&&( (rp2<=rd)||(abs(rp2-rd)/rd<=tol) )
	    push!(edge_dielectric, k)
	    push!(edge_dielec_ind, findfirst(x->x==k, edgekeylist))
	end
    end 

    return edge_dielectric, edge_dielec_ind
end
