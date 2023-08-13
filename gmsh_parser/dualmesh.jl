# This file implements the functionality to create the complete dual mesh.


include("primalmesh.jl")


########################################### START DUAL NODES ###########################################
""" 
    get_circumcenter_tet(nodedict::Dict{Int, Nodestruct}, 
                         tet::Tetstruct)

Return the coords of the circumcenter of the tet.
"""
function get_circumcenter_tet(nodedict::Dict{Int, Nodestruct}, 
                              tet::Tetstruct)::SVector{3, Float64}

    nodeids = tet.nodes
    n1 = nodedict[nodeids[1]].coords
    n2 = nodedict[nodeids[2]].coords
    n3 = nodedict[nodeids[3]].coords
    n4 = nodedict[nodeids[4]].coords

    A = [transpose(n2-n1);
         transpose(n3-n1);
         transpose(n4-n1)]

    B = 0.5*[norm(n2)^2 - norm(n1)^2;
             norm(n3)^2 - norm(n1)^2;
             norm(n4)^2 - norm(n1)^2]

    return inv(A)*B
    
end


"""
    create_interior_dualnodedict(nodedict::Dict{Int, Nodestruct}, 
                                 tetdict::Dict{Int, Tetstruct})

Create the interior_dualnodedict, given tetdict.
"""
function create_interior_dualnodedict(nodedict::Dict{Int, Nodestruct}, 
                                      tetdict::Dict{Int, Tetstruct})::Dict{Int, Interior_dualnodestruct}

    interior_dualnodedict = Dict{Int, Interior_dualnodestruct}()

    for tetpair in tetdict

        tetid = tetpair.first
        tetstruct = tetpair.second

        # make interior_dualnodestruct & insert into dict
        interior_dualnode = Interior_dualnodestruct()
        interior_dualnode.id = tetid
        interior_dualnode.coords = get_circumcenter_tet(nodedict, tetstruct)
        
        interior_dualnodedict[tetid] = interior_dualnode
    
    end

    return interior_dualnodedict

end


"""
    get_tetsforface(face::Facestruct, 
                    tetdict::Dict{Int, Tetstruct})

Return the tetids that belong to face.
- Length 1 vector if face is on boundary
- length 2 if face is in interior (asc order of nodeid)
"""
function get_tetsforface(face::Facestruct, 
                         tetdict::Dict{Int, Tetstruct})::Vector{Int}

    # get set of nodes that make up face
    faceid = face.id

    # find which tets contain this face
    parent_tetids = []

    for tetpair in tetdict
        
        tetnodes = tetpair.second.nodes
        if issubset(faceid, tetnodes)
            push!(parent_tetids, tetpair.first)
        end
    
    end
    
    return sort(parent_tetids)

end


"""
    get_circumcenter_face(face::Facestruct, 
                               nodedict::Dict{Int, Nodestruct})

Return the circumcenter of face.
"""
function get_circumcenter_face(face::Facestruct, 
                               nodedict::Dict{Int, Nodestruct})::SVector{3, Float64}
    
    # get nodeids of face
    nodeids = face.id

    # get node coords
    a = nodedict[nodeids[1]].coords
    b = nodedict[nodeids[2]].coords
    c = nodedict[nodeids[3]].coords

    # get circumcenter
    circumcenter =  begin
                        a + 
                        ((norm(c-a)^2 * cross(cross(b-a, c-a), b-a)) 
                        + (norm(b-a))^2 * cross(cross(c-a, b-a), c-a)) /
                        (2 * norm(cross(b-a, c-a))^2)
                    end 
    
                    return circumcenter 

end


"""
    create_boundary_dualnodedict(nodedict::Dict{Int, Nodestruct},
                                      facedict::Dict{SVector{3, Int}, Facestruct},
                                      tetdict::Dict{Int, Tetstruct})

Create the boundary_dualnodedict, given facedict.
A face is on the boundary iff it belongs to only 1 tet.
"""
function create_boundary_dualnodedict(nodedict::Dict{Int, Nodestruct},
                                      facedict::Dict{SVector{3, Int}, Facestruct},
                                      tetdict::Dict{Int, Tetstruct})::Dict{SVector{3, Int}, Boundary_dualnodestruct}

    boundary_dualnodedict = Dict{SVector{3, Int}, Boundary_dualnodestruct}()

    for facepair in facedict

        faceid = facepair.first
        facestruct = facepair.second

        # get tets sharing face
        tetids = get_tetsforface(facestruct, tetdict)

        # if only one tet shares face, then face on boundary
        if length(tetids) == 1

            # make boundary dual node & insert into dict
            boundary_dualnode = Boundary_dualnodestruct()
            boundary_dualnode.id = faceid
            boundary_dualnode.coords = get_circumcenter_face(facestruct, nodedict)
            boundary_dualnode.tet = tetids[1]

            boundary_dualnodedict[faceid] = boundary_dualnode
        
        end
    
    end

    return boundary_dualnodedict

end


"""
    get_boundaryfaces_foredge(edge::Edgestruct, 
                                   boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})

Get the boundary faces that contain edge.
boundary_dualnodedict from previous function contains the dict of all possible boundary faces ids.

If length(get_boundaryfaces_foredge()) >= 1, then the edge belongs on the boundary.
"""
function get_boundaryfaces_foredge(edge::Edgestruct, 
                                   boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})::Vector{SVector{3, Int}}

    # get set of nodes that make up edge
    edgeid = edge.id

    # find which boundary faces contain this edge
    parent_boundary_faceids = []

    for boundary_faceid in keys(boundary_dualnodedict)

        if issubset(edgeid, boundary_faceid)

            push!(parent_boundary_faceids, boundary_faceid)

        end
    
    end
    
    return parent_boundary_faceids

end


"""
    gget_midpoint_edge(edge::Edgestruct, 
                           nodedict::Dict{Int, Nodestruct})

Return the midpoint of edge.
"""
function get_midpoint_edge(edge::Edgestruct, 
                           nodedict::Dict{Int, Nodestruct})::SVector{3, Float64}

    # get nodes
    edgeid = edge.id

    # return midpoint
    node1coords = nodedict[edgeid[1]].coords
    node2coords = nodedict[edgeid[2]].coords

    return 0.5 * (node1coords + node2coords)

end


"""
    create_auxiliary_dualnodedict(nodedict::Dict{Int, Nodestruct}, 
                                       edgedict::Dict{SVector{2, Int}, Edgestruct},
                                       boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})

Create the Auxiliary_dualnodedict, given edgedict.
"""
function create_auxiliary_dualnodedict(nodedict::Dict{Int, Nodestruct}, 
                                       edgedict::Dict{SVector{2, Int}, Edgestruct},
                                       boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})::Dict{SVector{2, Int}, Auxiliary_dualnodestruct}

    auxiliary_dualnodedict = Dict{SVector{2, Int}, Auxiliary_dualnodestruct}()

    for edgepair in edgedict

        edgeid = edgepair.first
        edgestruct = edgepair.second

        # get boundary faces sharing edge
        boundary_faceids = get_boundaryfaces_foredge(edgestruct, boundary_dualnodedict)

        # If length(get_boundaryfaces_foredge(edge)) >= 1, then the edge belongs on the boundary
        if length(boundary_faceids) >= 1

            # make auxiliary dual node & insert into dict
            auxiliary_dualnode = Auxiliary_dualnodestruct()
            auxiliary_dualnode.id = edgeid
            auxiliary_dualnode.coords = get_midpoint_edge(edgestruct, nodedict)

            auxiliary_dualnodedict[edgeid] = auxiliary_dualnode
        
        end
    
    end

    return auxiliary_dualnodedict

end


"""
    create_dualnodedicts(nodedict::Dict{Int, Nodestruct}, 
                              edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                              facedict::Dict{SVector{3, Int}, Facestruct}, 
                              tetdict::Dict{Int, Tetstruct})

Create dualnodedicts.
"""
function create_dualnodedicts(nodedict::Dict{Int, Nodestruct}, 
                              edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                              facedict::Dict{SVector{3, Int}, Facestruct}, 
                              tetdict::Dict{Int, Tetstruct})::Dualnodedicts_struct

    dualnodedicts = Dualnodedicts_struct()

    dualnodedicts.interior_dualnodedict = create_interior_dualnodedict(nodedict, tetdict)
    dualnodedicts.boundary_dualnodedict = create_boundary_dualnodedict(nodedict, facedict, tetdict)
    dualnodedicts.auxiliary_dualnodedict = create_auxiliary_dualnodedict(nodedict, edgedict, dualnodedicts.boundary_dualnodedict)

    return dualnodedicts

end
########################################### END DUAL NODES ###########################################


########################################### START DUAL EDGES ###########################################
"""
    create_interior_dualedgedict(facedict::Dict{SVector{3, Int}, Facestruct},
                                 tetdict::Dict{Int, Tetstruct}, 
                                 boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})

Create dictionary of interior dual edges.

The list of interior primal faces is obtained by the set difference between facedict and boundary_dualnodedict.
"""
function create_interior_dualedgedict(facedict::Dict{SVector{3, Int}, Facestruct},
                                      tetdict::Dict{Int, Tetstruct}, 
                                      interior_dualnodedict::Dict{Int, Interior_dualnodestruct}, 
                                      boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})::Dict{SVector{3, Int}, Interior_dualedgestruct}

    interior_dualedgedict = Dict{SVector{3, Int}, Interior_dualedgestruct}()

    # get list of interior primal faces 
    interior_primalfaceids = setdiff(keys(facedict), keys(boundary_dualnodedict))
    
    for interior_primalfaceid in interior_primalfaceids

        interior_dualedge = Interior_dualedgestruct()
        interior_dualedge.id = interior_primalfaceid
        interior_dualedge.dualnodes = get_tetsforface(facedict[interior_primalfaceid], tetdict) # [interior dual node id 1, interior dual node id 2], in ascending order

        vec = interior_dualnodedict[interior_dualedge.dualnodes[1]].coords - interior_dualnodedict[interior_dualedge.dualnodes[2]].coords
        interior_dualedge.length = norm(vec)

        # insert into dict
        interior_dualedgedict[interior_dualedge.id] = interior_dualedge

    end

    return interior_dualedgedict

end



"""
    create_boundary_dualedgedict(interior_dualnodedict::Dict{Int, Interior_dualnodestruct}, 
                                 boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})

Create dictionary of boundary dual edges. 

The list of boundary primal face ids are contained in boundary_dualnodedict.
"""
function create_boundary_dualedgedict(interior_dualnodedict::Dict{Int, Interior_dualnodestruct}, 
                                      boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct})::Dict{SVector{3, Int}, Boundary_dualedgestruct}

    boundary_dualedgedict = Dict{SVector{3, Int}, Boundary_dualedgestruct}()                                  

    for boundary_dualnodepair in boundary_dualnodedict

        boundary_dualnodeid = boundary_dualnodepair.first
        boundary_dualnode = boundary_dualnodepair.second

        boundary_dualedge = Boundary_dualedgestruct()
        boundary_dualedge.id = boundary_dualnodeid
        boundary_dualedge.dualnodes = [boundary_dualnode.tet, boundary_dualnodeid]  # [interior dual node id, boundary dual node id] 
        
        vec = interior_dualnodedict[boundary_dualedge.dualnodes[1]].coords - boundary_dualnodedict[boundary_dualedge.dualnodes[2]].coords
        boundary_dualedge.length = norm(vec)

        # insert into dict
        boundary_dualedgedict[boundary_dualedge.id] = boundary_dualedge
         
    end

    return boundary_dualedgedict

end


"""
    create_auxiliary_onprimaledge_dualedgedict(nodedict::Dict{Int, Nodestruct},
                                               auxiliary_dualnodedict)

Create dictionary of auxiliary dual edges lying on part of a boundary primal edge.

Each such dual edge corresponds to a tuple (boundary primal edge, primal node part of that boundary primal edge).
Boundary primal edge ids are already listed in auxiliary_dualnodedict.
"""
function create_auxiliary_onprimaledge_dualedgedict(nodedict::Dict{Int, Nodestruct},
                                                    auxiliary_dualnodedict)::Dict{SVector{2, Any}, Auxiliary_onprimaledge_dualedgestruct}

    auxiliary_onprimaledge_dualedgedict = Dict{SVector{2, Any}, Auxiliary_onprimaledge_dualedgestruct}()           
    
    # make auxiliary_onprimaledge_dualedgedict, by looping over each tuple (boundary primal edge, primal node part of that boundary primal edge), 
    # & insert into dict
    for boundary_primaledgeid in keys(auxiliary_dualnodedict)

        for boundary_primalnodeid in boundary_primaledgeid

            auxiliary_onprimaledge_dualedge = Auxiliary_onprimaledge_dualedgestruct()
            auxiliary_onprimaledge_dualedge.id = [boundary_primaledgeid, boundary_primalnodeid] 
            
            vec = auxiliary_dualnodedict[auxiliary_onprimaledge_dualedge.id[1]].coords - nodedict[auxiliary_onprimaledge_dualedge.id[2]].coords
            auxiliary_onprimaledge_dualedge.length = norm(vec)

            auxiliary_onprimaledge_dualedgedict[auxiliary_onprimaledge_dualedge.id] = auxiliary_onprimaledge_dualedge

        end

    end

    return auxiliary_onprimaledge_dualedgedict

end



"""
    create_auxiliary_onprimalface_dualedgedict(facedict::Dict{SVector{3, Int}, Facestruct},
                                               boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct}, 
                                               auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct})

Create dictionary of auxiliary dual edges lying on the area of a boundary primal face.

Each such dual edge corresponds to a tuple (boundary primal face, primal edge part of that boundary primal face).
Boundary primal face ids are already listed in boundary_dualnodedict.
"""
function create_auxiliary_onprimalface_dualedgedict(facedict::Dict{SVector{3, Int}, Facestruct},
                                                    boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct},
                                                    auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct})::Dict{SVector{2, Any}, Auxiliary_onprimalface_dualedgestruct}

    auxiliary_onprimalface_dualedgedict = Dict{SVector{2, Any}, Auxiliary_onprimalface_dualedgestruct}()           
    
    # make auxiliary_onprimalface_dualedgedict, by looping over each tuple (boundary primal face, primal edge part of that boundary primal face), 
    # & insert into dict
    for boundary_primalfaceid in keys(boundary_dualnodedict)

        for boundary_primaledgeid in facedict[boundary_primalfaceid].edges

            auxiliary_onprimalface_dualedge = Auxiliary_onprimalface_dualedgestruct()
            auxiliary_onprimalface_dualedge.id = [boundary_primalfaceid, boundary_primaledgeid] 
            
            vec = boundary_dualnodedict[auxiliary_onprimalface_dualedge.id[1]].coords - auxiliary_dualnodedict[auxiliary_onprimalface_dualedge.id[2]].coords
            auxiliary_onprimalface_dualedge.length = norm(vec)

            auxiliary_onprimalface_dualedgedict[auxiliary_onprimalface_dualedge.id] = auxiliary_onprimalface_dualedge

        end

    end

    return auxiliary_onprimalface_dualedgedict

end


"""
    create_dualedgedicts(nodedict::Dict{Int, Nodestruct}, 
                         facedict::Dict{SVector{3, Int}, Facestruct}, 
                         tetdict::Dict{Int, Tetstruct}, 
                         interior_dualnodedict::Dict{Int, Interior_dualnodestruct}, 
                         boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct},
                         auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct})


Create dualedgedicts.
"""
function create_dualedgedicts(nodedict::Dict{Int, Nodestruct}, 
                              facedict::Dict{SVector{3, Int}, Facestruct}, 
                              tetdict::Dict{Int, Tetstruct}, 
                              interior_dualnodedict::Dict{Int, Interior_dualnodestruct}, 
                              boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct},
                              auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct})
    
    dualedgedicts = Dualedgedicts_struct()
    dualedgedicts.interior_dualedgedict = create_interior_dualedgedict(facedict, tetdict, interior_dualnodedict, boundary_dualnodedict)
    dualedgedicts.boundary_dualedgedict = create_boundary_dualedgedict(interior_dualnodedict, boundary_dualnodedict)
    dualedgedicts.auxiliary_onprimaledge_dualedgedict = create_auxiliary_onprimaledge_dualedgedict(nodedict, auxiliary_dualnodedict)
    dualedgedicts.auxiliary_onprimalface_dualedgedict = create_auxiliary_onprimalface_dualedgedict(facedict, boundary_dualnodedict, auxiliary_dualnodedict)

    return dualedgedicts

end
########################################### END DUAL EDGES ###########################################


########################################### START DUAL FACES ###########################################
"""
    get_tetsforedge_unordered(edge::Edgestruct, 
                              tetdict::Dict{Int, Tetstruct})

Return the vector of tets which contain edge.
Note however that this vector is not ordered properly.
"""
function get_tetsforedge_unordered(edge::Edgestruct, 
                                   tetdict::Dict{Int, Tetstruct})::Vector{Int}

    # get set of nodes that make up edge
    edgeid = edge.id

    # find which tets contain this edge
    parent_tetids = []

    for tetpair in tetdict
        tetnodes = tetpair.second.nodes
        if issubset(edgeid, tetnodes)
            push!(parent_tetids, tetpair.first)
        end

    end

    return parent_tetids

end


"""
    get_nexttet(current_tetid::Int, 
                         last_tetid::Int, 
                         parent_tetids::Vector{Int}, 
                         tetdict::Dict{Int, Tetstruct})

Return the next tetid starting at current_tetid, 
determined using parent_tetids (all tet_ids shared by edge)
and excluding last_tetid.
"""
function get_nexttet(current_tetid::Int, 
                     last_tetid::Int, 
                     parent_tetids::Vector{Int}, 
                     tetdict::Dict{Int, Tetstruct})

    adjacent_tetids = []

    current_tetnodes = tetdict[current_tetid].nodes

    # see which tet_ids in parent_tetids are adjacent to current tet
    for parent_tetid in parent_tetids

        adjacent = false 

        # want to exclude current_tetid from consideration
        if parent_tetid != current_tetid
            
            parent_tetnodes = tetdict[parent_tetid].nodes
            
            # if shared face, append to list of adjacent tetids
            if length(intersect(current_tetnodes, parent_tetnodes)) == 3

                push!(adjacent_tetids, parent_tetid)
            
            end

        end

    end

    # the next tet is the one with min id,
    # and want to exclude last_tetid from consideration

    next_tet_vec = setdiff(adjacent_tetids, last_tetid)
    if length(next_tet_vec) != 0
        return next_tet_vec[1]
    else
        return Set() # return empty set if no such tets exist
    end

end

# when possible, RENAME THIS get_tetsforprimaledge_ordered
"""
get_tetsforedge_ordered(edge::Edgestruct, 
                        tetdict::Dict{Int, Tetstruct})


Return the vector of tets which contain (an interior primal) edge, ordered.
"""
function get_tetsforedge_ordered(edge::Edgestruct, 
                                 tetdict::Dict{Int, Tetstruct})::Vector{Int}

    parent_tetids_ordered = []

    # list of all tets containing this edge
    parent_tetids_unordered = get_tetsforedge_unordered(edge, tetdict::Dict{Int, Tetstruct})

    # start with lowest tet id
    low_tetid = minimum(parent_tetids_unordered)
    push!(parent_tetids_ordered, low_tetid)

    # call get_nexttet() once to begin the sequence
    current_tetid = get_nexttet(low_tetid, low_tetid, parent_tetids_unordered, tetdict)
    next_tetid = current_tetid
    last_tetid = low_tetid

    # get the next tet
    while (next_tetid != low_tetid) && length(next_tetid) != 0

        push!(parent_tetids_ordered, next_tetid)
        next_tetid = get_nexttet(current_tetid, last_tetid, parent_tetids_unordered, tetdict)
        last_tetid = current_tetid
        current_tetid = next_tetid

    end

    return parent_tetids_ordered

end


"""
get_tetsforboundaryedge_ordered(edge::Edgestruct, 
                        tetdict::Dict{Int, Tetstruct})


Return the vector of tets which contain a boundary edge, ordered.
"""
function get_tetsforboundaryedge_ordered(edge::Edgestruct, 
                                         tetdict::Dict{Int, Tetstruct},
                                         boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})::Vector{Int}

    parent_tetids_ordered = []

    # list of all tets containing this edge
    parent_tetids_unordered = get_tetsforedge_unordered(edge, tetdict::Dict{Int, Tetstruct})

    # determine the subset of 2 tets which contain the 2 boundary faces
    tets_containing_boundary_primalfaces = []
    
    boundary_primalfaces = get_boundaryfaces_for_boundaryedge(edge, boundary_dualedgedict)
    
    for parent_tetid in parent_tetids_unordered
        faces = tetdict[parent_tetid].faces
        for face in faces
            # check if any of the tet's faces match the 2 boundary faces
            if face in boundary_primalfaces
                push!(tets_containing_boundary_primalfaces, parent_tetid)
            end
        end
    end

    # start with lowest tet id OUT OF the two tets which contain the 2 boundary primal faces
    low_tetid = minimum(tets_containing_boundary_primalfaces)
    push!(parent_tetids_ordered, low_tetid)

    # call get_nexttet() once to begin the sequence
    current_tetid = get_nexttet(low_tetid, low_tetid, parent_tetids_unordered, tetdict)
    next_tetid = current_tetid
    last_tetid = low_tetid

    # get the next tet
    while (next_tetid != low_tetid) && length(next_tetid) != 0

        push!(parent_tetids_ordered, next_tetid)
        next_tetid = get_nexttet(current_tetid, last_tetid, parent_tetids_unordered, tetdict)
        last_tetid = current_tetid
        current_tetid = next_tetid

    end

    return parent_tetids_ordered

end


"""

Get the interior dual edges belonging to interior dual face.
"""
function get_interior_dualedges_for_interior_dualface(interior_dualface::Interior_dualfacestruct,
                                                      interior_dualedgedict::Dict{SVector{3, Int}, Interior_dualedgestruct})

    interior_dualnodes = interior_dualface.interior_dualnodes

    # want list of interior dual edge ids (i.e. interior primal face ids)
    desired_interior_dualedgeids = []

    for i in 1:lastindex(interior_dualnodes)-1

        # interior_dualedgedict lists dualnodes for each interior dual edge in ascending order
        interior_dualedge_candidate = sort([interior_dualnodes[i], interior_dualnodes[i+1]])

        for interior_dualedgepair in interior_dualedgedict

            interior_dualedgestruct = interior_dualedgepair.second

            interior_dualedge_dualnodes = interior_dualedgestruct.dualnodes

            if interior_dualedge_candidate == interior_dualedge_dualnodes
                push!(desired_interior_dualedgeids, interior_dualedgepair.first)
                break
            end

        end

    end 

    # last interior dual edge needs to wrap around the end to the start of interior_dualnodes
    interior_dualedge_candidate = sort([interior_dualnodes[end], interior_dualnodes[1]])
    
    for interior_dualedgepair in interior_dualedgedict

        interior_dualedgestruct = interior_dualedgepair.second

        interior_dualedge_dualnodes = interior_dualedgestruct.dualnodes

        if interior_dualedge_candidate == interior_dualedge_dualnodes
            push!(desired_interior_dualedgeids, interior_dualedgepair.first)
            break
        end

    end

    return desired_interior_dualedgeids

end



"""


Return the unsigned triangle area made of a 2D elementary dual entity.

We use the word "triangle" to simply distinguish the elementary 
dual entity from the triangles in the primal mesh.
"""
function get_area_triangle(coords1::SVector{3, Float64}, 
                           coords2::SVector{3, Float64}, 
                           coords3::SVector{3, Float64})::Float64

    area = 0.5*abs(norm(cross(coords2-coords1, coords3-coords1)))

end


"""

Determine the raw value of the dual area using Hirani's method.
"""
function get_dualarea_rawvalue(edge::Edgestruct,
                               nodedict::Dict{Int, Nodestruct},
                               facedict::Dict{SVector{3, Int}, Facestruct}, 
                               tetdict::Dict{Int, Tetstruct})::Float64

    dualarea_value = 0

    # get the list of faces containing edge
    edgeid = edge.id

    facelist = []
    for facepair in facedict
        faceid = facepair.first
        if issubset(edgeid, faceid)
            push!(facelist, faceid)
        end
    end
    
    # for each face in above list,
    # get list of tets containing that face
    for faceid in facelist

        face = facedict[faceid]

        tetlist = []
        for tetpair in tetdict 
            tetid = tetpair.first
            tetnodes = tetpair.second.nodes
            if issubset(faceid, tetnodes)
                push!(tetlist, tetid)
            end
        end

        for tetid in tetlist

            tet = tetdict[tetid]

            # calculate multiplier

            # calculate lambda1 (omitted, since it is always 1)
            mp = get_midpoint_edge(edge, nodedict)

            lambda1 = 1

            # calculate lambda2 
            vp2 = setdiff(faceid,  edgeid)[1]

            vp2_coords = nodedict[vp2].coords

            vector1 = vp2_coords - nodedict[edgeid[1]].coords
            vector2 = vp2_coords - nodedict[edgeid[2]].coords
            face_normal = cross(vector1, vector2)
            vector3 = nodedict[edgeid[2]].coords - nodedict[edgeid[1]].coords
            halfspace_plane_normal = cross(face_normal, vector3)

            a = halfspace_plane_normal[1]
            b = halfspace_plane_normal[2]
            c = halfspace_plane_normal[3]
            x0 = nodedict[edgeid[1]].coords[1]
            y0 = nodedict[edgeid[1]].coords[2]
            z0 = nodedict[edgeid[1]].coords[3]
            face_circumcenter = get_circumcenter_face(face, nodedict)

            sign_vp2 = sign(a*(vp2_coords[1] - x0) + b*(vp2_coords[2] - y0) + c*(vp2_coords[3] - z0))
            sign_face_circumcenter = sign(a*(face_circumcenter[1] - x0) + b*(face_circumcenter[2] - y0) + c*(face_circumcenter[3] - z0))

            lambda2 = -1
            if sign_vp2 == sign_face_circumcenter
                lambda2 = 1
            end

            # calculate lambda3
            vp3 = setdiff(tet.nodes, faceid)[1]
            tet_circumcenter = get_circumcenter_tet(nodedict, tetdict[tetid])

            vp3_coords = nodedict[vp3].coords

            d = face_normal[1]
            e = face_normal[2]
            f = face_normal[3]

            sign_vp3 = sign(d*(vp3_coords[1] - x0) + e*(vp3_coords[2] - y0) + f*(vp3_coords[3] - z0))
            sign_tet_circumcenter = sign(d*(tet_circumcenter[1] - x0) + e*(tet_circumcenter[2] - y0) + f*(tet_circumcenter[3] - z0))

            lambda3 = -1
            if sign_vp3 == sign_tet_circumcenter
                lambda3 = 1
            end

            # calculate overall multiplier
            multiplier = lambda1 * lambda2 * lambda3

            # add to dual volume
            triangle_area = get_area_triangle(mp, face_circumcenter, tet_circumcenter)
            dualarea_value += multiplier * triangle_area

        end

    end 

    return dualarea_value

end


"""

Determine the effective value of the dual area using Hirani's method.
"""
function get_dualarea_effectivevalue(edge::Edgestruct,
                               nodedict::Dict{Int, Nodestruct},
                               facedict::Dict{SVector{3, Int}, Facestruct}, 
                               tetdict::Dict{Int, Tetstruct}, material::Array{Float64,1})::Float64

    dualarea_value = 0

    # get the list of faces containing edge
    edgeid = edge.id

    facelist = []
    for facepair in facedict
        faceid = facepair.first
        if issubset(edgeid, faceid)
            push!(facelist, faceid)
        end
    end
    
    # for each face in above list,
    # get list of tets containing that face
    for faceid in facelist

        face = facedict[faceid]

        tetlist = []
        for tetpair in tetdict 
            tetid = tetpair.first
            tetnodes = tetpair.second.nodes
            if issubset(faceid, tetnodes)
                push!(tetlist, tetid)
            end
        end

        for tetid in tetlist

            tet = tetdict[tetid]

            # calculate multiplier

            # calculate lambda1 (omitted, since it is always 1)
            mp = get_midpoint_edge(edge, nodedict)

            lambda1 = 1

            # calculate lambda2 
            vp2 = setdiff(faceid,  edgeid)[1]

            vp2_coords = nodedict[vp2].coords

            vector1 = vp2_coords - nodedict[edgeid[1]].coords
            vector2 = vp2_coords - nodedict[edgeid[2]].coords
            face_normal = cross(vector1, vector2)
            vector3 = nodedict[edgeid[2]].coords - nodedict[edgeid[1]].coords
            halfspace_plane_normal = cross(face_normal, vector3)

            a = halfspace_plane_normal[1]
            b = halfspace_plane_normal[2]
            c = halfspace_plane_normal[3]
            x0 = nodedict[edgeid[1]].coords[1]
            y0 = nodedict[edgeid[1]].coords[2]
            z0 = nodedict[edgeid[1]].coords[3]
            face_circumcenter = get_circumcenter_face(face, nodedict)

            sign_vp2 = sign(a*(vp2_coords[1] - x0) + b*(vp2_coords[2] - y0) + c*(vp2_coords[3] - z0))
            sign_face_circumcenter = sign(a*(face_circumcenter[1] - x0) + b*(face_circumcenter[2] - y0) + c*(face_circumcenter[3] - z0))

            lambda2 = -1
            if sign_vp2 == sign_face_circumcenter
                lambda2 = 1
            end

            # calculate lambda3
            vp3 = setdiff(tet.nodes, faceid)[1]
            tet_circumcenter = get_circumcenter_tet(nodedict, tetdict[tetid])

            vp3_coords = nodedict[vp3].coords

            d = face_normal[1]
            e = face_normal[2]
            f = face_normal[3]

            sign_vp3 = sign(d*(vp3_coords[1] - x0) + e*(vp3_coords[2] - y0) + f*(vp3_coords[3] - z0))
            sign_tet_circumcenter = sign(d*(tet_circumcenter[1] - x0) + e*(tet_circumcenter[2] - y0) + f*(tet_circumcenter[3] - z0))

            lambda3 = -1
            if sign_vp3 == sign_tet_circumcenter
                lambda3 = 1
            end

            # calculate overall multiplier
            multiplier = lambda1 * lambda2 * lambda3

            # add to dual volume
            triangle_area = get_area_triangle(mp, face_circumcenter, tet_circumcenter)
            dualarea_value += multiplier * triangle_area * material[tetdict[tetid].entityid]

        end

    end 

    return dualarea_value

end


"""
    get_suportvolume(edge::Edgestruct, dualfacedicts::Dualfacedicts_struct)::Float64

Compute the raw support volumes of primal edges
"""
function get_suportvolume(edge::Edgestruct, dualfacedicts::Dualfacedicts_struct)::Float64

    interioredge_key = keys(dualfacedicts.interior_dualfacedict)
    boundaryedge_key = keys(dualfacedicts.boundary_dualfacedict) 

    edge_id = edge.id
    if edge_id in interioredge_key
        supportvolume_raw = dualfacedicts.interior_dualfacedict[edge_id].raw_area*edge.length
    elseif edge_id in boundaryedge_key
        supportvolume_raw = dualfacedicts.boundary_dualfacedict[edge_id].raw_area*edge.length
    end

    return supportvolume_raw
    
end


"""

Create dictionary of interior dual faces.

Boundary primal edge ids are already listed in auxiliary_dualnodedict,
so we take the set setdiff(keys(edgedict), keys(auxiliary_dualnodedict)) to find the interior primal edge ids
"""
function create_interior_dualfacedict(nodedict::Dict{Int, Nodestruct},
                                      edgedict::Dict{SVector{2, Int}, Edgestruct},
                                      facedict::Dict{SVector{3, Int}, Facestruct},
                                      tetdict::Dict{Int, Tetstruct},
                                      auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct},
                                      interior_dualedgedict::Dict{SVector{3, Int}, Interior_dualedgestruct})::Dict{SVector{2, Int}, Interior_dualfacestruct}

    interior_dualfacedict = Dict{SVector{2, Int}, Interior_dualfacestruct}()

    # Boundary primal edge ids are already listed in auxiliary_dualnodedict,
    # so we take the set setdiff(keys(edgedict), keys(auxiliary_dualnodedict)) to find the interior primal edge ids
    for interior_edgeid in setdiff(keys(edgedict), keys(auxiliary_dualnodedict))

        interior_dualface = Interior_dualfacestruct()
        interior_dualface.id = interior_edgeid
        interior_dualface.interior_dualnodes = get_tetsforedge_ordered(edgedict[interior_edgeid], tetdict)
        interior_dualface.interior_dualedges = get_interior_dualedges_for_interior_dualface(interior_dualface, interior_dualedgedict)
        interior_dualface.raw_area = get_dualarea_rawvalue(edgedict[interior_edgeid], nodedict, facedict, tetdict)

        # insert into dict
        interior_dualfacedict[interior_dualface.id] = interior_dualface

    end

    return interior_dualfacedict

end



"""
get_boundaryfaces_for_boundaryedge(edge::Edgestruct, 
                                   boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})

Return the vector of boundary faces which contain boundary_edge, ordered lexicographically. 
This should return a length 2 vector.

The list of boundary faces is in boundary_dualedgedict.
"""
function get_boundaryfaces_for_boundaryedge(edge::Edgestruct, 
                                            boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})::Vector{SVector{3, Int}}

    edgeid = edge.id

    # find which faces contain this edge
    parent_faceids = []

    for facepair in boundary_dualedgedict
        faceid = facepair.first
        if issubset(edgeid, faceid)
            push!(parent_faceids, faceid)
        end
    end

    # sort lexicographically
    return sort(parent_faceids)

end


"""

Get the interior dual edges belonging to boundary dual face.
"""
function get_interior_dualedges_for_boundary_dualface(boundary_dualface::Boundary_dualfacestruct,
                                                      interior_dualedgedict::Dict{SVector{3, Int}, Interior_dualedgestruct})::Vector{SVector{3, Int}}

    # boundary_dual face has dualnodes: [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]  
    # only want the interior ones                                             
    interior_dualnodes = boundary_dualface.dualnodes[3:end-1] 

    # want list of interior dual edge ids (i.e. interior primal face ids)
    desired_interior_dualedgeids = []

    # want list of interior dual edge ids (i.e. interior primal face ids)
    desired_interior_dualedgeids = []

    for i in 1:lastindex(interior_dualnodes)-1

        # interior_dualedgedict lists dualnodes for each interior dual edge in ascending order
        interior_dualedge_candidate = sort([interior_dualnodes[i], interior_dualnodes[i+1]])

        for interior_dualedgepair in interior_dualedgedict

            interior_dualedgestruct = interior_dualedgepair.second

            interior_dualedge_dualnodes = interior_dualedgestruct.dualnodes

            if interior_dualedge_candidate == interior_dualedge_dualnodes
                push!(desired_interior_dualedgeids, interior_dualedgepair.first)
                break
            end

        end

    end 

    # NO WRAPPING AROUND 
    # last interior dual edge needs to wrap around the end to the start of interior_dualnodes
    # interior_dualedge_candidate = sort!([interior_dualnodes[end], interior_dualnodes[1]])
    # 
    # for interior_dualedgepair in interior_dualedgedict
    #
    #    interior_dualedgestruct = interior_dualedgepair.second
    #
    #    interior_dualedge_dualnodes = interior_dualedgestruct.dualnodes
    #
    #        if interior_dualedge_candidate == interior_dualedge_dualnodes
    #            push!(desired_interior_dualedgeids, interior_dualedgepair.first)
    #            break
    #        end
    #
    # end

    return desired_interior_dualedgeids

end


"""

Get the boundary dual edges belonging to boundary dual face.
"""
function get_boundary_dualedges_for_boundary_dualface(boundary_dualface::Boundary_dualfacestruct, 
                                                      boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})::Vector{SVector{3, Int}}

    # dual nodes = [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]                                                  
    dualnodes = boundary_dualface.dualnodes

    # want list of boundary dual edge ids (i.e. boundary primal face ids)
    desired_boundary_dualedgeids = []

    for boundary_dualedgepair in boundary_dualedgedict

        boundary_dualedge_id = boundary_dualedgepair.first
        boundary_dualedge = boundary_dualedgepair.second

        boundary_dualedge_dualnodes = boundary_dualedge.dualnodes

        if boundary_dualedge_dualnodes == [dualnodes[3], dualnodes[2]] # [interior dual node 1, boundary dual node 1] 
            push!(desired_boundary_dualedgeids, boundary_dualedge_id)
            break
        end

    end

    
    for boundary_dualedgepair in boundary_dualedgedict

        boundary_dualedge_id = boundary_dualedgepair.first
        boundary_dualedge = boundary_dualedgepair.second

        boundary_dualedge_dualnodes = boundary_dualedge.dualnodes

        if boundary_dualedge_dualnodes == [dualnodes[end-1], dualnodes[end]] # [interior dual node n-1, boundary dual node 2] 
            push!(desired_boundary_dualedgeids, boundary_dualedge_id)
            break
        end

    end

    return desired_boundary_dualedgeids

end


"""

Get the auxiliary_onprimalface_dualedges belonging to boundary dual face.
"""
function get_auxiliary_onprimalface_dualedges_for_boundary_dualface(boundary_dualface::Boundary_dualfacestruct)::Vector{SVector{2, Any}}

    # dual nodes = [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]                                                  
    dualnodes = boundary_dualface.dualnodes

    # want list of auxiliary_onprimalface_dualedges ids (i.e. boundary primal edge ids)
    desired_auxiliary_onprimalface_dualedge_ids = []

    push!(desired_auxiliary_onprimalface_dualedge_ids, [dualnodes[2], dualnodes[1]]) # [boundary dual node id 1, auxiliary dual node id] 
    push!(desired_auxiliary_onprimalface_dualedge_ids, [dualnodes[end], dualnodes[1]]) # [boundary dual node id 2, auxiliary dual node id] 

    return desired_auxiliary_onprimalface_dualedge_ids

end


"""

Create dictionary of boundary dual faces.

Boundary primal edge ids are already listed in auxiliary_dualnodedict.
"""
function create_boundary_dualfacedict(nodedict::Dict{Int, Nodestruct},
                                      edgedict::Dict{SVector{2, Int}, Edgestruct},
                                      facedict::Dict{SVector{3, Int}, Facestruct},
                                      tetdict::Dict{Int, Tetstruct},
                                      auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct},
                                      interior_dualedgedict::Dict{SVector{3, Int}, Interior_dualedgestruct},
                                      boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})::Dict{SVector{2, Int}, Boundary_dualfacestruct}

    boundary_dualfacedict = Dict{SVector{2, Int}, Boundary_dualfacestruct}()

    # Boundary primal edge ids are already listed in auxiliary_dualnodedict
    for boundary_edgeid in keys(auxiliary_dualnodedict)

        boundary_dualface =  Boundary_dualfacestruct()
        boundary_dualface.id = boundary_edgeid

        boundaryfaces = get_boundaryfaces_for_boundaryedge(edgedict[boundary_edgeid], boundary_dualedgedict)

        # get list of tets containing edge,
        # starting with the tet containing face 1 
        # and ending with the tet containing face 2
        # (if needed, flip the default order) #
        parent_tetids = get_tetsforboundaryedge_ordered(edgedict[boundary_edgeid], tetdict, boundary_dualedgedict) 
        if parent_tetids[1] != get_tetsforface(facedict[boundaryfaces[1]], tetdict)[1]
            reverse!(parent_tetids)
        end

        # [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]
        dualnodes = parent_tetids
        dualnodes = convert(Vector{Any}, dualnodes) # convert to Vector{Any} since we immediately inserting objects other than Ints
        insert!(dualnodes, 1, boundary_edgeid)
        insert!(dualnodes, 2, boundaryfaces[1])
        push!(dualnodes, boundaryfaces[2])

        boundary_dualface.dualnodes = dualnodes

        boundary_dualface.interior_dualedges = get_interior_dualedges_for_boundary_dualface(boundary_dualface, interior_dualedgedict)

        boundary_dualface.boundary_dualedges = get_boundary_dualedges_for_boundary_dualface(boundary_dualface, boundary_dualedgedict)
        
        boundary_dualface.auxiliary_onprimalface_dualedges = get_auxiliary_onprimalface_dualedges_for_boundary_dualface(boundary_dualface)
        
        boundary_dualface.raw_area = get_dualarea_rawvalue(edgedict[boundary_edgeid], nodedict, facedict, tetdict)

        # insert into dict
        boundary_dualfacedict[boundary_dualface.id] = boundary_dualface
        
    end

    return boundary_dualfacedict

end


"""

Create dictionary of auxiliary dual faces.

Corresponds bijectively to a tuple (boundary primal face, primal node part of that primal face).
The boundary faces are contained in boundary_dualedgedict.
"""
function create_auxiliary_dualfacedict(facedict::Dict{SVector{3, Int}, Facestruct},
                                       boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct})::Dict{SVector{2, Any} , Auxiliary_dualfacestruct}

    auxiliary_dualfacedict = Dict{SVector{2, Any} , Auxiliary_dualfacestruct}()

    for boundary_faceid in keys(boundary_dualedgedict)

        for boundary_primalnodeid in boundary_faceid

            auxiliary_dualface = Auxiliary_dualfacestruct()
            
            auxiliary_dualface.id = [boundary_faceid, boundary_primalnodeid] # [boundary primal face id, primal node part of that primal face]

            # get edges of face containing node
            parent_edges = []
            face_edges = facedict[boundary_faceid].edges
            for face_edge in face_edges
                if boundary_primalnodeid in face_edge
                    push!(parent_edges, face_edge)
                end
            end
            sort!(parent_edges)
            
            auxiliary_dualface.dualnodes = [boundary_primalnodeid, parent_edges[1], boundary_faceid, parent_edges[2]] #[primal node, auxiliary dual node 1, boundary dual node, auxiliary dual node 2]
            

            # get auxiliary_onprimaledge_dualedges for auxiliary_dualface
            # each such dual edge is of the form [boundary primal edge id, primal node part of that boundary primal edge id]
            auxiliary_onprimaledge_dualedges = []
            push!(auxiliary_onprimaledge_dualedges, [parent_edges[1], boundary_primalnodeid])
            push!(auxiliary_onprimaledge_dualedges, [parent_edges[2], boundary_primalnodeid])

            auxiliary_dualface.auxiliary_onprimaledge_dualedges = auxiliary_onprimaledge_dualedges 

            # get auxiliary_onprimalface_dualedges for auxiliary_dualface
            # each such dual edge is of the form [boundary primal face id, boundary primal edge id]
            auxiliary_onprimalface_dualedges = []
            push!(auxiliary_onprimalface_dualedges, [boundary_faceid, parent_edges[1]]) 
            push!(auxiliary_onprimalface_dualedges, [boundary_faceid, parent_edges[2]])

            auxiliary_dualface.auxiliary_onprimalface_dualedges = auxiliary_onprimalface_dualedges 

            # insert into dict
            auxiliary_dualfacedict[auxiliary_dualface.id] = auxiliary_dualface

        end

    end

    return auxiliary_dualfacedict

end
########################################### END DUAL FACES ###########################################


########################################### START DUAL VOLUMES ###########################################
"""

Return the unsigned pyramid volume.

A "pyramid" refers to the tetrahedron formed
by the coordinates of a single 3D elementary dual entity.
Its nodes are the circumcenters of the following:
[node_id, edge_id, face_id, tet_id].

We use the word "pyramid" simply to distinguish the elementary dual entity
from the tetrahedra in the primal mesh.
"""
function get_volume_pyramid(coords1::SVector{3, Float64}, 
                            coords2::SVector{3, Float64}, 
                            coords3::SVector{3, Float64}, 
                            coords4::SVector{3, Float64})::Float64

    volume = 1/6*abs(dot(coords4-coords1, cross(coords2-coords1, coords3-coords1)))

end


"""

Determine the raw value of the dual volume using Hirani's method.
"""
function get_dualvolume_rawvalue(node::Nodestruct, 
                                 nodedict,
                                 edgedict, 
                                 facedict, 
                                 tetdict)::Float64

    dualvolume_value = 0

    nodeid = node.id

    # get list of edges containing node
    edgelist = []
    for edgepair in edgedict
        edgeid = edgepair.first
        if nodeid in edgeid
            push!(edgelist, edgeid)
        end
    end

    # for each edge in above list, 
    # get the list of faces containing that edge
    for edgeid in edgelist

        edge = edgedict[edgeid]

        facelist = []
        for facepair in facedict
            faceid = facepair.first
            if issubset(edgeid, faceid)
                push!(facelist, faceid)
            end
        end
        
        # for each face in above list,
        # get list of tets containing that face
        for faceid in facelist

            face = facedict[faceid]

            tetlist = []
            for tetpair in tetdict 
                tetid = tetpair.first
                tetnodes = tetpair.second.nodes
                if issubset(faceid, tetnodes)
                    push!(tetlist, tetid)
                end
            end

            for tetid in tetlist

                tet = tetdict[tetid]

                # calculate multiplier

                # calculate lambda1 (omitted, since it is always 1)
                mp = get_midpoint_edge(edge, nodedict)

                lambda1 = 1

                # calculate lambda2 
                vp2 = setdiff(faceid,  edgeid)[1]

                vp2_coords = nodedict[vp2].coords

                vector1 = vp2_coords - nodedict[edgeid[1]].coords
                vector2 = vp2_coords - nodedict[edgeid[2]].coords
                face_normal = cross(vector1, vector2)
                vector3 = nodedict[edgeid[2]].coords - nodedict[edgeid[1]].coords
                halfspace_plane_normal = cross(face_normal, vector3)

                a = halfspace_plane_normal[1]
                b = halfspace_plane_normal[2]
                c = halfspace_plane_normal[3]
                x0 = nodedict[edgeid[1]].coords[1]
                y0 = nodedict[edgeid[1]].coords[2]
                z0 = nodedict[edgeid[1]].coords[3]
                face_circumcenter = get_circumcenter_face(face, nodedict)

                sign_vp2 = sign(a*(vp2_coords[1] - x0) + b*(vp2_coords[2] - y0) + c*(vp2_coords[3] - z0))
                sign_face_circumcenter = sign(a*(face_circumcenter[1] - x0) + b*(face_circumcenter[2] - y0) + c*(face_circumcenter[3] - z0))

                lambda2 = -1
                if sign_vp2 == sign_face_circumcenter
                    lambda2 = 1
                end

                # calculate lambda3
                vp3 = setdiff(tet.nodes, faceid)[1]
                tet_circumcenter = get_circumcenter_tet(nodedict, tetdict[tetid])

                vp3_coords = nodedict[vp3].coords

                d = face_normal[1]
                e = face_normal[2]
                f = face_normal[3]

                sign_vp3 = sign(d*(vp3_coords[1] - x0) + e*(vp3_coords[2] - y0) + f*(vp3_coords[3] - z0))
                sign_tet_circumcenter = sign(d*(tet_circumcenter[1] - x0) + e*(tet_circumcenter[2] - y0) + f*(tet_circumcenter[3] - z0))

                lambda3 = -1
                if sign_vp3 == sign_tet_circumcenter
                    lambda3 = 1
                end

                # calculate overall multiplier
                multiplier = lambda1 * lambda2 * lambda3

                # add to dual volume
                pyramid_volume = get_volume_pyramid(node.coords, mp, face_circumcenter, tet_circumcenter)
                dualvolume_value += multiplier * pyramid_volume

            end

        end 

    end

    return dualvolume_value

end


"""


Return the vector of edges which contain node.
This will be used in get_dualvolume() to get the list of dual faces corresponding to primal node.
"""
function get_edgesfornode(node::Nodestruct, 
                          edgedict::Dict{SVector{2, Int}, Edgestruct})::Vector{SVector{2, Int}}

    nodeid = node.id

    # find which edges contain this tet
    parent_edgeids = []

    for edgepair in edgedict
        edgeid = edgepair.first
        if issubset(nodeid, edgeid)
            push!(parent_edgeids, edgeid)
        end

    end

    return parent_edgeids

end


"""


Return the dual volume corresponding to node.
"""
function get_dualvolume(node::Nodestruct, 
                        nodedict::Dict{Int, Nodestruct},
                        edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                        facedict::Dict{SVector{3, Int}, Facestruct}, 
                        tetdict::Dict{Int, Tetstruct},
                        interior_dualfacedict::Dict{SVector{2, Int}, Interior_dualfacestruct},  
                        boundary_dualfacedict::Dict{SVector{2, Int}, Boundary_dualfacestruct}, 
                        auxiliary_dualfacedict::Dict{SVector{2, Any}, Auxiliary_dualfacestruct})::Dualvolumestruct
    
    dualvolume = Dualvolumestruct()
    
    dualvolume.raw_volume = get_dualvolume_rawvalue(node, nodedict, edgedict, facedict, tetdict)
    
    dualvolume.interior_dualfaces = Vector{SVector{2, Int}}()
    dualvolume.boundary_dualfaces = Vector{SVector{2, Int}}()
    dualvolume.auxiliary_dualfaces = Vector{SVector{2, Any}}()

    # get edges containing node
    parent_edgeids = get_edgesfornode(node, edgedict)

    # see which dict each edge above belongs to
    for parent_edgeid in parent_edgeids

        if parent_edgeid in keys(boundary_dualfacedict)

            push!(dualvolume.boundary_dualfaces, parent_edgeid)
        
        else

            push!(dualvolume.interior_dualfaces, parent_edgeid)
        
        end
    
    end

    # append auxiliary dual face, if its id [boundary primal face id, primal node part of that primal face] contains node
    for auxiliary_dualface_pair in auxiliary_dualfacedict

        auxiliary_dualface_id = auxiliary_dualface_pair.first
        nodeid = node.id

        if nodeid == auxiliary_dualface_id[2]

            push!(dualvolume.auxiliary_dualfaces, auxiliary_dualface_id)

        end 

    end

    # complete rest of info for dual volume

    # get all interior dual nodes, by looping over each of the dual volume's
    # interior dual face
    # and retrieving the interior dual nodes
    # use union of sets to only include unique elements 
    # (NOTE: the advantage of this method, say compared to determining all tets containing node, 
    # is that we don't need to loop over all tets).
    all_interior_dualnodes = Set([])
    for interior_dualfaceid in dualvolume.interior_dualfaces   
        union!(all_interior_dualnodes, Set(interior_dualfacedict[interior_dualfaceid].interior_dualnodes))
    end

    dualvolume.interior_dualnodes = collect(all_interior_dualnodes) #convert back to vector

    # get all boundary dual nodes, by looping over each of the dual volume's
    # boundary_dualface [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]
    # and each auxiliary_dualface [primal node, auxiliary dual node 1, boundary dual node, auxiliary dual node 2],
    # and retrieving the boundary dual nodes
    all_boundary_dualnodes = Set([])
    
    for boundary_dualfaceid in dualvolume.boundary_dualfaces
        union!(all_boundary_dualnodes, Set(boundary_dualfacedict[boundary_dualfaceid].dualnodes[[2, end]]))
    end
    for auxiliary_dualfaceid in dualvolume.auxiliary_dualfaces
        union!(all_boundary_dualnodes, Set(auxiliary_dualfacedict[auxiliary_dualfaceid].dualnodes[[3]]))
    end
    

    # println(all_boundary_dualnodes)
    dualvolume.boundary_dualnodes = collect(all_boundary_dualnodes)
    
    # get auxiliary dual nodes, by looping over each of the dual volume's
    # boundary_dualface [auxiliary dual node, boundary dual node 1, interior dual nodes ..., boundary dual node 2]
    # and each auxiliary_dualface [primal node, auxiliary dual node 1, boundary dual node, auxiliary dual node 2],
    # and retrieving the auxiliary dual nodes
    all_auxiliary_dualnodes = Set([])

    for boundary_dualfaceid in dualvolume.boundary_dualfaces
        union!(all_auxiliary_dualnodes, Set(boundary_dualfacedict[boundary_dualfaceid].dualnodes[[1]]))
    end
    
    for auxiliary_dualfaceid in dualvolume.auxiliary_dualfaces
        union!(all_auxiliary_dualnodes, Set(auxiliary_dualfacedict[auxiliary_dualfaceid].dualnodes[[2, 4]]))
    end
    

    dualvolume.auxiliary_dualnodes = collect(all_auxiliary_dualnodes)

    # get interior_dualedges, by looping over each of the dual volume's
    # interior dual face
    # and each boundary_dualface
    # and retrieving the interior dual edges
    all_interior_dualedges = Set([])
    for interior_dualfaceid in dualvolume.interior_dualfaces 
        union!(all_interior_dualedges, Set(interior_dualfacedict[interior_dualfaceid].interior_dualedges))
    end
    for boundary_dualfaceid in dualvolume.boundary_dualfaces 
        union!(all_interior_dualedges, Set(boundary_dualfacedict[boundary_dualfaceid].interior_dualedges))
    end 

    dualvolume.interior_dualedges = collect(all_interior_dualedges)

    # get boundary_dualedges, by looping over each of the dual volume's
    # boundary_dualface
    # and retrieving the boundary dual edges
    all_boundary_dualedges = Set([])
    for boundary_dualfaceid in dualvolume.boundary_dualfaces 
        union!(all_boundary_dualedges, Set(boundary_dualfacedict[boundary_dualfaceid].boundary_dualedges))
    end

    dualvolume.boundary_dualedges = collect(all_boundary_dualedges)

    # get auxiliary_onprimalface_dualedges, by looping over each of the dual volume's
    # boundary_dualface
    # and auxiliary_dualface
    # and retrieving the auxiliary_onprimalface_dualedges
    all_auxiliary_onprimalface_dualedges = Set([])
    for  boundary_dualfaceid in dualvolume.boundary_dualfaces 
        union!(all_auxiliary_onprimalface_dualedges, Set(boundary_dualfacedict[boundary_dualfaceid].auxiliary_onprimalface_dualedges))
    end
    for  auxiliary_dualfaceid in dualvolume.auxiliary_dualfaces 
        union!(all_auxiliary_onprimalface_dualedges, Set(auxiliary_dualfacedict[auxiliary_dualfaceid].auxiliary_onprimalface_dualedges))
    end

    dualvolume.auxiliary_onprimalface_dualedges = collect(all_auxiliary_onprimalface_dualedges)

    # get auxiliary_onprimaledge_dualedges, by looping over each of the dual volume's
    # auxiliary_dualface
    # and retrieving the auxiliary_onprimaledge_dualedges
    all_auxiliary_onprimaledge_dualedges = Set([])
    for  auxiliary_dualfaceid in dualvolume.auxiliary_dualfaces 
        union!(all_auxiliary_onprimaledge_dualedges, Set(auxiliary_dualfacedict[auxiliary_dualfaceid].auxiliary_onprimaledge_dualedges))
    end

    dualvolume.auxiliary_onprimaledge_dualedges = collect(all_auxiliary_onprimaledge_dualedges)

    return dualvolume

end


"""
create_dualvolumedict(nodedict::Dict{Int, Nodestruct}, 
                      edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                      facedict::Dict{SVector{3, Int}, Facestruct}, 
                      tetdict::Dict{Int, Tetstruct},
                      interior_dualfacedict::Dict{SVector{2, Int}, Interior_dualfacestruct},  
                      boundary_dualfacedict::Dict{SVector{2, Int}, Boundary_dualfacestruct}, 
                      auxiliary_dualfacedict::Dict{SVector{2, Any}, Auxiliary_dualfacestruct})

Create dictionary of dual volumes.
"""
function create_dualvolumedict(nodedict::Dict{Int, Nodestruct}, 
                               edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                               facedict::Dict{SVector{3, Int}, Facestruct}, 
                               tetdict::Dict{Int, Tetstruct},
                               interior_dualfacedict::Dict{SVector{2, Int}, Interior_dualfacestruct},  
                               boundary_dualfacedict::Dict{SVector{2, Int}, Boundary_dualfacestruct}, 
                               auxiliary_dualfacedict::Dict{SVector{2, Any}, Auxiliary_dualfacestruct})::Dict{Int, Dualvolumestruct}

    dualvolumedict = Dict{Int, Dualvolumestruct}()

    for nodepair in nodedict

        nodeid = nodepair.first
        node = nodepair.second

        dualvolume = get_dualvolume(node, 
                                    nodedict,
                                    edgedict, 
                                    facedict, 
                                    tetdict, 
                                    interior_dualfacedict,
                                    boundary_dualfacedict,
                                    auxiliary_dualfacedict)

        dualvolumedict[nodeid] = dualvolume

    end

    return dualvolumedict

end
########################################### END DUAL VOLUMES ###########################################


########################################### START DUAL MESH ###########################################
"""
    complete_dualmesh(file::String)

Returns all completed dual mesh dictionaries.
"""
function complete_dualmesh(file::String)

    primalmesh, _, _ = complete_primalmesh(file)

    nodedict = primalmesh.nodedict
    edgedict = primalmesh.edgedict
    facedict = primalmesh.facedict
    tetdict = primalmesh.tetdict

    dualmesh = Dualmeshstruct()

    # create dualnodedicts & insert into dual mesh
    dualnodedicts = Dualnodedicts_struct()
    dualnodedicts.interior_dualnodedict = create_interior_dualnodedict(nodedict, 
                                                                       tetdict)
    dualnodedicts.boundary_dualnodedict = create_boundary_dualnodedict(nodedict, 
                                                                       facedict, 
                                                                       tetdict)
    dualnodedicts.auxiliary_dualnodedict = create_auxiliary_dualnodedict(nodedict, 
                                                                         edgedict, 
                                                                         dualnodedicts.boundary_dualnodedict)

    dualmesh.dualnodedicts =  dualnodedicts

    # create dualedgedicts & insert into dual mesh
    dualedgedicts = Dualedgedicts_struct()
    dualedgedicts.interior_dualedgedict = create_interior_dualedgedict(facedict, 
                                                                       tetdict,
                                                                       dualnodedicts.interior_dualnodedict,
                                                                       dualnodedicts.boundary_dualnodedict)
    dualedgedicts.boundary_dualedgedict = create_boundary_dualedgedict(dualnodedicts.interior_dualnodedict, 
                                                                       dualnodedicts.boundary_dualnodedict)

    dualedgedicts.auxiliary_onprimaledge_dualedgedict = create_auxiliary_onprimaledge_dualedgedict(nodedict, 
                                                                                                   dualnodedicts.auxiliary_dualnodedict)
    dualedgedicts.auxiliary_onprimalface_dualedgedict = create_auxiliary_onprimalface_dualedgedict(facedict, 
                                                                                                   dualnodedicts.boundary_dualnodedict, 
                                                                                                   dualnodedicts.auxiliary_dualnodedict)

    dualmesh.dualedgedicts = dualedgedicts

    # create dualfacedicts & insert into dual mesh
    dualfacedicts = Dualfacedicts_struct()
    dualfacedicts.interior_dualfacedict = create_interior_dualfacedict(nodedict,
                                                                       edgedict,
                                                                       facedict,
                                                                       tetdict,
                                                                       dualnodedicts.auxiliary_dualnodedict,
                                                                       dualedgedicts.interior_dualedgedict)
    dualfacedicts.boundary_dualfacedict = create_boundary_dualfacedict(nodedict,
                                                                       edgedict,
                                                                       facedict,
                                                                       tetdict,
                                                                       dualnodedicts.auxiliary_dualnodedict,
                                                                       dualedgedicts.interior_dualedgedict,
                                                                       dualedgedicts.boundary_dualedgedict)
    dualfacedicts.auxiliary_dualfacedict = create_auxiliary_dualfacedict(facedict, 
                                                                         dualedgedicts.boundary_dualedgedict)

    dualmesh.dualfacedicts = dualfacedicts

    # create dualvomedict & insert into dual mesh
    dualmesh.dualvolumedict = create_dualvolumedict(nodedict,
                                                    edgedict,
                                                    facedict,
                                                    tetdict,
                                                    dualfacedicts.interior_dualfacedict,
                                                    dualfacedicts.boundary_dualfacedict,
                                                    dualfacedicts.auxiliary_dualfacedict)                       


    # compute raw support volumes for primal edges
    for edgeid in keys(edgedict)
        edgedict[edgeid].supportvolume = get_suportvolume(edgedict[edgeid], dualfacedicts)
    end
    # update edgedict
    primalmesh.edgedict = edgedict

    return dualmesh, primalmesh

end
########################################### END DUAL MESH ###########################################


#= note, on Windows, use \ instead of /

@time primalmesh = complete_primalmesh(raw"/meshes/transfinite_test.msh")
nodedict = primalmesh.nodedict
edgedict = primalmesh.edgedict
facedict = primalmesh.facedict
tetdict = primalmesh.tetdict


dualnodedicts = Dualnodedicts_struct()
dualnodedicts.interior_dualnodedict = create_interior_dualnodedict(nodedict, 
                                                                    tetdict)
dualnodedicts.boundary_dualnodedict = create_boundary_dualnodedict(nodedict, 
                                                                    facedict, 
                                                                    tetdict)
dualnodedicts.auxiliary_dualnodedict = create_auxiliary_dualnodedict(nodedict, 
                                                                     edgedict, 
                                                                     dualnodedicts.boundary_dualnodedict)

dualedgedicts = Dualedgedicts_struct()
dualedgedicts.interior_dualedgedict = create_interior_dualedgedict(facedict, 
                                                                    tetdict,
                                                                    dualnodedicts.interior_dualnodedict,
                                                                    dualnodedicts.boundary_dualnodedict)
dualedgedicts.boundary_dualedgedict = create_boundary_dualedgedict(dualnodedicts.interior_dualnodedict, 
                                                                    dualnodedicts.boundary_dualnodedict)

# [199, [6, 14, 43]] works

# [174, [1, 17, 32]] for some reason, is not appearing as dualnodes field in boundary_dualedgedict

# PROBLEM: none of the faces of tet 174 is are at the boundary, so [174, [1, 17, 32]] should not be a candidate in first place
# Any[[1, 32], [1, 9, 32], 205, 174, [1, 17, 32]]
# should be 
# Any[[1, 32], [1, 9, 32], 205, 193, [1, 17, 32]]

# CAUSE: get_tetsforedge_ordered() is buggy
# this is because starting at the min tet id out of ALL tet ids sharing edge
# DOES NOT WORK for boundary dual edges, as the loop will not be complete;
# currently only works for interior dual edges.
# Need a different convention for tet id ordering for boundary dual edges, for example,
# the first tet id is the min out of the 2 tets corresponding to the 2 boundary faces.

#println(keys(dualnodedicts.boundary_dualnodedict))
"""

# @time dualmesh = complete_dualmesh(raw"/meshes/transfinite_test.msh")

=#