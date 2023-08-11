function pp_sharp(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, edgeform::Array{Float64,1}, Nv::Int64, Ne::Int64)
    
    pvec = Array{Float64,2}(undef, Nv, 3)
    nodedict = primalmesh.nodedict
    edgedict = primalmesh.edgedict
    facedict = primalmesh.facedict
    tetdict  = primalmesh.tetdict
    nodekeylist = collect(keys(nodedict))
    edegkeylist = collect(keys(edgedict))
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
            edge_ind = findfirst(x->x==edgeid, edegkeylist)
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


function pp_sharp_lq(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, edgeform::Array{Float64,1}, Nv::Int64, Ne::Int64)

end
