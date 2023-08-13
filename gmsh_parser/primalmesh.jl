# This file implements the functionality to create the complete primal mesh.

include("parse.jl")
using LinearAlgebra
using Combinatorics


"""
    gettetvol(tet:Tetstruct)

Return the volume of tet.

Uses the cross-product formula.
"""
function gettetvol(nodedict::Dict{Int, Nodestruct}, tet::Tetstruct)::Float64

    # get nodeids belonging to tet
    nodeids = tet.nodes

    # get coords of nodes
    node1_coords = nodedict[nodeids[1]].coords
    node2_coords = nodedict[nodeids[2]].coords
    node3_coords = nodedict[nodeids[3]].coords
    node4_coords = nodedict[nodeids[4]].coords

    # use cross-product formula to get volume
    a = node2_coords - node1_coords
    b = node3_coords - node1_coords
    c = node4_coords - node1_coords

    # return volume
    volume = 1/6 * abs(dot(cross(a, b), c))

end


"""
    getedgelen(edge::Edgestruct)

Return the length of edge.
"""
function getedgelen(nodedict::Dict{Int, Nodestruct}, edge::Edgestruct)::Float64
    
    # get coords of nodes
    node1coords = nodedict[edge.id[1]].coords
    node2coords = nodedict[edge.id[2]].coords

    # calculate length
    vec = node1coords - node2coords
    length = norm(vec)
    return length
end


"""
    getfacearea(face::Facestruct)

Return the area of face.
"""
function getfacearea(nodedict::Dict{Int, Nodestruct}, face::Facestruct)

    # get coords of nodes
    node1coords = nodedict[face.id[1]].coords
    node2coords = nodedict[face.id[2]].coords
    node3coords = nodedict[face.id[3]].coords

    # calculate area using cross product
    vec1 = node2coords - node1coords
    vec2 = node3coords - node1coords
    area = 1/2 * norm(cross(vec1, vec2))
    return area
end


"""
    complete_update_tet(tet::Tetstruct, edgedict, facedict)

Updates the missing fields of tet and updates the associated structs
to the appropriate primal mesh dictionaries.
"""
function complete_update_tet(nodedict::Dict{Int, Nodestruct}, edgedict::Dict{SVector{2, Int}, Edgestruct}, facedict::Dict{SVector{3, Int}, Facestruct}, tet::Tetstruct)

    # get nodeids
    nodeids = tet.nodes

    # instantiate empty Vectors for tet
    edgeid_tet_vec = Vector{}()
    faceid_tet_vec = Vector{}()

    faceid_tet_vec = sort!.(collect(combinations(nodeids, 3)))

    for faceid in faceid_tet_vec
        
        # get 3 edgeids
        edgeids = sort!.(collect(combinations(faceid, 2)))

        for edgeid in edgeids
            
            # make Edgestruct & get length, then append to edgedict
            edge = Edgestruct()
            edge.id = edgeid
            edge.length = getedgelen(nodedict, edge)
            #edge.entities_dict[1] = [0]
            #edge.entities_dict[2] = [0]
            #edge.entities_dict[3] = [0]
            edgedict[edgeid] = edge

            # append edgeid to tet, if not already in
            if ~(edgeid in edgeid_tet_vec)
                push!(edgeid_tet_vec, edgeid)
            end

        end

        # make Facestruct & get area, then append to facedict
        face = Facestruct()
        face.id = faceid
        face.edges = edgeids
        face.area = getfacearea(nodedict, face)
        facedict[faceid] = face
        
    end
    
    tet.edges = edgeid_tet_vec
    tet.faces = faceid_tet_vec
    tet.volume = gettetvol(nodedict, tet)

end


"""
    edge2entities_map(nodedict::Dict{Int, Nodestruct}, edgedict::Dict{SVector{2, Int}, Edgestruct})

adding information about which mesh entities each edge belongs to
"""
function edge2entities_map(nodedict::Dict{Int, Nodestruct}, edgedict::Dict{SVector{2, Int}, Edgestruct})

    for ekey in keys(edgedict)
        n1 = ekey[1]
        n2 = ekey[2]
        entities_dict = Dict{Int,Vector{}}()
        if (nodedict[n1].entities_dict[1]!=[0])&&(nodedict[n2].entities_dict[1]!=[0])
            # Note: and edge can only lie on either 1 or 0 zero physical curves, hence findfirst
            sharedcurve_ind = findfirst(in(nodedict[n1].entities_dict[1]), nodedict[n2].entities_dict[1])
            if !isnothing(sharedcurve_ind) # if n1 and n2 both lie on the same curve
                entities_dict[1] = [nodedict[n2].entities_dict[1][sharedcurve_ind]]   
            else
                entities_dict[1] = [0]           
            end
        else 
            entities_dict[1] = [0]           
        end

        if (nodedict[n1].entities_dict[2]!=[0])&&(nodedict[n2].entities_dict[2]!=[0])
            # Note: if an edge lies on a curve, that edge can belong to multiple neighboring surfaces, hence findall
            sharedsurface_ind = findall(in(nodedict[n1].entities_dict[2]), nodedict[n2].entities_dict[2])
            if !isnothing(sharedsurface_ind) # if n1 and n2 both lie on the same surface
                entities_dict[2] = nodedict[n2].entities_dict[2][sharedsurface_ind]
            else
                entities_dict[2] = [0]
            end
        else
            entities_dict[2] = [0]
        end

        if (nodedict[n1].entities_dict[3]!=[0])&&(nodedict[n2].entities_dict[3]!=[0])
            # Note: if an edge lies on a curve or a surface, that edge can belong to multiple neighboring volumes, hence findall
            sharedvolume_ind = findall(in(nodedict[n1].entities_dict[3]), nodedict[n2].entities_dict[3])
            if !isnothing(sharedvolume_ind) # if n1 and n2 both lie on the same volume
                entities_dict[3] = nodedict[n2].entities_dict[3][sharedvolume_ind]
            else
                entities_dict[3] = [0]
            end
        else
            entities_dict[3] = [0]
        end

        edgedict[ekey].entities_dict = entities_dict
    end
    #return edgedict
end

"""
    complete_primalmesh(file::String)

Returns the 4 completed primal mesh dictionaries.
"""
function complete_primalmesh(file::String)


    #nodedict, tetdict = parsefile(file)
    nodedict, tetdict, physicalnames_dict, all_entities_struct = parsefile(file)
    edgedict = Dict{SVector{2, Int}, Edgestruct}()
    facedict = Dict{SVector{3, Int}, Facestruct}()
    
    for tetpair in tetdict
        complete_update_tet(nodedict, edgedict, facedict, tetpair.second)
    end

    edge2entities_map(nodedict, edgedict)

    primalmesh = Primalmeshstruct()
    primalmesh.nodedict = nodedict
    primalmesh.edgedict = edgedict
    primalmesh.facedict = facedict
    primalmesh.tetdict  = tetdict

    return primalmesh, physicalnames_dict, all_entities_struct
end


# @time primalmesh = complete_primalmesh(raw"\meshes\transfinite_test.msh")
# 
# println("number of nodes: ", length(primalmesh[1]))
# println("number of edges: ", length(primalmesh[2]))
# println("number of faces: ", length(primalmesh[3]))
# println("number of tets: ", length(primalmesh[4]))

