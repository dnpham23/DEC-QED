# compute the primal-primal sharp operator
function pp_sharp(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, edgeform::Array{Float64,1}, Nv::Int64, Ne::Int64)
    
    pvec = Array{Float64,2}(undef, Nv, 3)
    nodedict = primalmesh.nodedict
    edgedict = primalmesh.edgedict
    facedict = primalmesh.facedict
    tetdict  = primalmesh.tetdict
    nodekeylist = collect(keys(nodedict))
    edgekeylist = collect(keys(edgedict))
    tetkeylist  = collect(keys(tetdict))
    
    for i = 1:Nv
        nodeid  = nodekeylist[i]
        tempvec = [0.0, 0.0, 0.0]
        
        # get list of edges containing node
        edgelist = []
        for edge in edgedict
            edgeid = edge.first
            if nodeid in edgeid
                push!(edgelist, edgeid)
            end
        end
        
        for edgeid in edgelist
            w = [0.0, 0.0, 0.0]
            # find the other node in this edge that's not the node in consideration
            node2 = edgeid[findfirst(x->x!=nodeid, edgeid)] 
            # find the tetrahedra that contain this edge
            for tetkey in tetkeylist
                if (edgeid in tetdict[tetkey].edges)
                    # find the face in the tet that's opposite to node2
                    opposite_face = tetdict[tetkey].faces[findfirst(x-> !(node2 in x), tetdict[tetkey].faces)]
                    # find the normal vector of this face
                    e1   = facedict[opposite_face].edges[1]
                    vec1 = nodedict[e1[1]].coords - nodedict[e1[2]].coords                                       
                    e2   = facedict[opposite_face].edges[2]
                    vec2 = nodedict[e2[1]].coords - nodedict[e2[2]].coords
                    
                    facenormal = cross(vec1,vec2)
                    facenormal_normalized = facenormal./norm(facenormal)
                    
                    # make sure the normal vector points to the same half-space that has node2
                    localvec = nodedict[node2].coords - nodedict[e1[1]].coords
                    facenormal_normalized = facenormal_normalized*sign(dot(facenormal_normalized, localvec))
                    
                    # find the gradient of the primal-primal interpolation function
                    h = 3.0*tetdict[tetkey].volume/facedict[opposite_face].area
                    gradInterp = 1/h*facenormal_normalized
                    
                    # find the three edges of the tet that contain the node in consideration
                    shared_edge_list = tetdict[tetkey].edges[findall(x->nodeid in x, tetdict[tetkey].edges)] 
                    
                    # compute the ratio of total dual volume inside the primal tet
                    dualvolume_intersect = 0.0
                    for s_edge in shared_edge_list
                        edge_circumcenter  = 0.5*(nodedict[s_edge[1]].coords + nodedict[s_edge[2]].coords)
                        # find the two faces of the tet that contain this edge
                        shared_face_list = tetdict[tetkey].faces[findall(x->(s_edge[1] in x)&&(s_edge[2] in x), tetdict[tetkey].faces)]
                        for s_face in shared_face_list
                            face_circumcenter     = get_circumcenter_face(facedict[s_face], nodedict)
                            tet_circumcenter      = get_circumcenter_tet(nodedict, tetdict[tetkey])
                            dualvolume_intersect += get_volume_pyramid(nodedict[nodeid].coords, edge_circumcenter, face_circumcenter, tet_circumcenter)
                        end
                    end
                    volume_ratio = dualvolume_intersect/tetdict[tetkey].volume
                    w           += volume_ratio*gradInterp
                end
            end
            
            # find the location of the current edge in the list of 1-form fields
            edge_ind = findfirst(x->x==edgeid, edgekeylist)
            # find the orientation of the current edge
            if nodeid == edgeid[1]
                sgn = 1.0
            else
                sgn = -1.0
            end
            
            # add the total contribution of this edge to the vector field in the node in consideration
            tempvec += sgn*edgeform[edge_ind]*w 
        end
        pvec[i,:] = tempvec
    end
    
    return pvec
end

# compute the primal-dual sharp operator to map the primal 1-forms on the boundary surface r3 to the tangential vector field that lives on the 
# circumcenters of the primal triangles on that surface
function pd_sharp_surface3(primalmesh::Primalmeshstruct, edgeform::Array{Float64,1}, ebound3::Array{Any,1}, fbound3::Array{Any,1})
    
    Ntri_bound3 = length(fbound3)
    dvec = Array{Float64,2}(undef, Ntri_bound3, 3)
    edgekeylist = collect(keys(edgedict))
    
    for nt = 1:Ntri_bound3
        fkey = fbound3[nt]
        edgelist = primalmesh.facedict[fkey].edges
        circ = get_circumcenter_face(facedict[fkey], nodedict)

        temp_dvec = zeros(3)
        for vert in fkey
            # find the other two vertices in the triangle
            v1 = fkey[findall(x->x!= vert, fkey)][1]
            v2 = fkey[findall(x->x!= vert, fkey)][2]
            # find the edges opposite to v1 and v2
            op_edge1 = edgelist[findfirst(x->!(v1 in x), edgelist)]
            op_edge2 = edgelist[findfirst(x->!(v2 in x), edgelist)]
            # height from vertices to opposite edges
            h1 = 2.0*primalmesh.facedict[fkey].area/primalmesh.edgedict[op_edge1].length
            h2 = 2.0*primalmesh.facedict[fkey].area/primalmesh.edgedict[op_edge2].length
            # find the vector normal to the opposite edge and points towards the vertex
            edge_center_1 = 0.5*(primalmesh.nodedict[op_edge1[1]].coords + primalmesh.nodedict[op_edge1[2]].coords) 
            edge_center_2 = 0.5*(primalmesh.nodedict[op_edge2[1]].coords + primalmesh.nodedict[op_edge2[2]].coords)
            edge2circ_1 = circ - edge_center_1 # assuming circ is inside the triangle
            edge2circ_2 = circ - edge_center_2 # assuming circ is inside the triangle
            e2c_len1     = norm(edge2circ_1)
            e2c_len2     = norm(edge2circ_2)
            n1 = edge2circ_1/c2e_len1 # normalized
            n2 = edge2circ_2/c2e_len2 # normalized
            # compute the weight based on the subdivision of the triangle
            w  = 0.25*(primalmesh.edgedict[op_edge1].length*e2c_len1 + primalmesh.edgedict[op_edge2].length*e2c_len2)/primalmesh.facedict[fkey].area
            # find the location of the edge [vert, v1] in the list of 1-form fields
            edge_ind1 = findfirst(x->x==op_edge2, edgekeylist) # note that the edge [vert, v1] is in fact op_edge2
            # find the orientation of the edge [vert, v1]
            if vert == op_edge2[1] # note that the edge [vert, v1] is in fact op_edge2
                sgn1 = 1.0
            else
                sgn1 = -1.0
            end
            # find the location of the edge [vert, v2] in the list of 1-form fields
            edge_ind2 = findfirst(x->x==op_edge1, edgekeylist) # note that the edge [vert, v2] is in fact op_edge1
            # find the orientation of the edge [vert, v2]
            if vert == op_edge1[1] # note that the edge [vert, v2] is in fact op_edge1
                sgn2 = 1.0
            else
                sgn2 = -1.0
            end
            temp_dvec += w*(sgn1*edgeform[edge_ind1]*(1/h1)*n1 + sgn2*edgeform[edge_ind2]*(1/h2)*n2)
        end
        dvec[nt] = temp_dvec
    end
    return dvec
end


# approximate the primal-primal sharp using least square method
function pp_sharp_least_sq_openbc(primalmesh::Primalmeshstruct, edgeform::Array{Float64,1}, edgekeylist::Array{SArray{Tuple{2},Int64,1,2},1}, 
                                  nodekeylist::Array{Any,1}, vbound4::Array{Any,1}, nodedict_extra::Dict{Int64,Nodestruct})
    Nv::Int64, Ne::Int64
    pvec = Array{Float64,2}(undef, Nv, 3) 

    
    for i = 1:Nv
        nodeid  = nodekeylist[i]
   
        # get list of edges containing node
        edgelist = []
        for edgeid in edgekeylist
            if nodeid == edgeid[1]
                push!(edgelist, [edgeid,1])
            elseif nodeid == edgeid[2]
                push!(edgelist, [edgeid,-1])
            end
        end

        nedge = length(edgelist)
        phi_v_list = zeros(nedge)
        edgevecs   = zeros(nedge,3)
        for ne = 1:nedge
            edge = edgelist[ne]
            edgeid = edge[1]
            v1 = edgeid[1]
            v2 = edgeid[2]
            if !(v1 in vbound4)
                v1_coord = primalmesh.nodedict[v1].coords
            else
                v1_coord = nodedict_extra[v1].coords
            end

            if !(v2 in vbound4)
                v2_coord = primalmesh.nodedict[v2].coords
            else
                v2_coord = nodedict_extra[v2].coords
            end
            # find the index of this edgeid in the list of all edges
            e_ind = findfirst(x->x==edgeid, edgekeylist)
            #push!(phi_v_list, edge[2]*edgeform[e_ind])
            phi_v_list[ne] = edge[2]*edgeform[e_ind]
            edgevecs[ne,:] = v2_coord-v1_coord
        end
        pvec[i,:] = edgevecs\phi_v_list
    end
    return pvec; 
end
