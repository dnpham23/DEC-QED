# This file defines the custom data types used throuhgout the program to represent objects in the primal and dual mesh.


# Using static arrays may improve performance for certain purposes
using StaticArrays

########################################### Start material information ###########################################
mutable struct Physicalname_struct
    physicaltag::Int
    dimension::Int
    name::String
    Physicalname_struct() = new()
end

# note, we have included only the relevant
# information of each entity for the purposes of DEC-QED
# a given entity may belong to multiple physical names
# which is why PhysicalTags::Vector{Int}
mutable struct Point_entity_struct
    point_tag::Int
    physicaltags::Vector{Int}
    Point_entity_struct() = new()
end 


mutable struct Curve_entity_struct
    curve_tag::Int
    physicaltags::Vector{Int}
    boundingpoints::Vector{Int}
    Curve_entity_struct() = new()
end 


mutable struct Surface_entity_struct
    surface_tag::Int
    physicaltags::Vector{Int}
    boundingcurves::Vector{Int}
    Surface_entity_struct() = new()
end 


mutable struct Volume_entity_struct
    volume_tag::Int
    physicaltags::Vector{Int}
    boundingsurfaces::Vector{Int}
    Volume_entity_struct() = new()
end 


mutable struct All_entities_struct
    point_entities_dict::Dict{Int, Point_entity_struct}
    curve_entities_dict::Dict{Int, Curve_entity_struct}
    surface_entities_dict::Dict{Int, Surface_entity_struct}
    volume_entities_dict::Dict{Int, Volume_entity_struct}
    All_entities_struct() = new()
end

########################################### End material information ###########################################

########################################### Start primal mesh ###########################################


# declare struct to contain node information
mutable struct Nodestruct
    id::Int 
    coords::SVector{3, Float64} 
    root_entitydim::Int # MshFileVersion 4.1
    root_entityid::Int
    entities_dict::Dict{Int,Vector{}}
    Nodestruct() = new()
end


# declare struct to contain edge information
mutable struct Edgestruct
    id::SVector{2, Int} # ordered by ascending node id 
    length::Float64
    entities_dict::Dict{Int,Vector{}}
    supportvolume::Float64 
    Edgestruct() = new()
end


# declare struct to contain face information
mutable struct Facestruct
    id::SVector{3, Int}
    edges::SVector{3, SVector{2, Int}} 
    area::Float64
    Facestruct() = new()
end


# declare struct to contain tet information
mutable struct Tetstruct
    id::Int
    nodes::SVector{4, Int} 
    edges::SVector{6, SVector{2, Int}} 
    faces::SVector{4, SVector{3, Int}} 
    entityid::Int
    volume::Float64
    Tetstruct() = new()
end


# declare struct to contain all primal mesh dicts
mutable struct Primalmeshstruct
    nodedict::Dict{Int, Nodestruct}
    edgedict::Dict{SVector{2, Int}, Edgestruct}
    facedict::Dict{SVector{3, Int}, Facestruct}
    tetdict::Dict{Int, Tetstruct}
    Primalmeshstruct() = new()
end


########################################### End primal mesh ###########################################
########################################### Start dual mesh ###########################################


# declare struct to contain interior dual node information
# an interior dual node corresponds to the circumcenter of a primal tet
mutable struct Interior_dualnodestruct
    id::Int # primal tet id
    coords::SVector{3, Float64} 
    Interior_dualnodestruct() = new()
end


# declare struct to contain boundary dual node information
# a boundary dual node is the circumcenter of a boundary primal face
# arising due to the truncation of the dual mesh at the primal mesh boundary
mutable struct Boundary_dualnodestruct
    id::SVector{3, Int} # boundary primal face id
    coords::SVector{3, Float64}
    tet::Int # primal tet id, this info is contained here mainly so it can be used in Boundary_dualedgestruct
    Boundary_dualnodestruct() = new()
end


# declare struct to contain auxiliary dual node information
# an auxiliary dual node is the midpoint (i.e. circumcenter) of a boundary primal edge
# arising due to the truncation of the dual mesh at the primal mesh boundary
mutable struct Auxiliary_dualnodestruct
    id::SVector{2, Int} # boundary primal edge id
    coords::SVector{3, Float64}
    Auxiliary_dualnodestruct() = new()
end


# declare struct to contain interior dual edge information
# an interior dual edge corresponds to an interior primal face
# (i.e. a primal face which is shared by 2 tets)
mutable struct Interior_dualedgestruct
    id::SVector{3, Int} # interior primal face id
    dualnodes::SVector{2, Int} # [interior dual node id 1, interior dual node id 2], in ascending order
    length::Float64
    Interior_dualedgestruct() = new()
end 


# declare struct to contain boundary dual edge information
# a boundary dual edge corresponds to a boundary primal face
# (i.e. a primal face which is shared by only 1 tet)
mutable struct Boundary_dualedgestruct
    id::SVector{3, Int} # boundary primal face id
    dualnodes::SVector{2, Any} # [interior dual node id, boundary dual node id] 
    length::Float64
    Boundary_dualedgestruct() = new()
end


# declare structs to contain auxiliary dual edge information
# there are two possible types of auxiliary edges
# each auxiliary edge can either: 
# 1) lie on part of boundary primal edge (in particular, is exactly half of boundary primal edge),
# such that it corresponds to a tuple (boundary primal edge, primal node part of that boundary primal edge), 
# so that each such boundary primal edge corresponds to 2 such auxiliary dual edges;
# OR 2) lie on the area of a boundary primal face,
# such that it corresponds to a tuple (boundary primal face, primal edge part of that boundary primal face),
# so that each such boundary primal face corresponds to 3 such auxiliary dual edges


# an auxiliary dual edge lying on a boundary primary edge
# corresponds to a tuple (boundary primal edge, primal node part of that boundary primal edge)
mutable struct Auxiliary_onprimaledge_dualedgestruct
    id::SVector{2, Any} # [boundary primal edge id, primal node part of that boundary primal edge id] , this is equivalent to its dualnodes [auxiliary dual node id, primal node part of that boundary primal edge id] 
    length::Float64
    Auxiliary_onprimaledge_dualedgestruct() = new()
end


# an auxiliary dual edge lying on a boundary primary face
# corresponds to a tuple (boundary primal face, primal edge part of that boundary primal face)
mutable struct Auxiliary_onprimalface_dualedgestruct
    id::SVector{2, Any} # [boundary primal face id, boundary primal edge id], this is equivalent to its dualnodes [boundary dual node id, auxiliary dual node id] 
    length::Float64
    Auxiliary_onprimalface_dualedgestruct() = new()
end


# declare struct to contain interior dual face information
# an interior dual face corresponds to an interior primal edge
# (i.e. a primal edge which does not belong to any boundary primal face)
mutable struct Interior_dualfacestruct
    id::SVector{2, Int} # interior primal edge id
    interior_dualnodes::Vector{Int} # list of interior dual node ids, starting with the min tet id and with the order defined by get_tetsforedge()
    interior_dualedges::Vector{SVector{3, Int}} # list of interior dual edge ids
    raw_area::Float64
    Interior_dualfacestruct() = new()
end


# declare struct to contain boundary dual face information
# a boundary dual face corresponds to a boundary primal edge
# (i.e. a primal edge which belongs to at least one boundary primal face;
# in particular, it will belong to exactly two primal boundary faces)
# note, these dual faces DO NOT have auxiliary_onprimaledge_dualedges
mutable struct Boundary_dualfacestruct
    id::SVector{2, Int} # primal edge id
    dualnodes::Vector{Any} # [auxiliary dual node id, boundary dual node id 1, interior dual node ids ..., boundary dual node id 2], boundary dual node id 1 is smaller than boundary dual node id 2 lexicographically
    interior_dualedges::Vector{SVector{3, Int}}
    boundary_dualedges::SVector{2, SVector{3, Int}}
    auxiliary_onprimalface_dualedges::SVector{2, SVector{2, Any}}
    raw_area::Float64
    Boundary_dualfacestruct() = new()
end


# declare struct to contain auxiliary dual face information
# an auxiliary dual face corresponds to a tuple (boundary primal face, primal node part of that primal face)
# the 2 auxiliary dual nodes are ordered lexicographically.
# note, these dual faces DO NOT have interior or boundary dual edges
# note, the areas for these elements are not included because as of now they are not necessary
mutable struct Auxiliary_dualfacestruct
    id::SVector{2, Any} # [boundary primal face id, primal node part of that primal face]
    dualnodes::Vector{Any} # [boundary primal node id, auxiliary dual node id 1, boundary dual node id, auxiliary dual node id 2], auxiliary dual node id 1 is smaller than auxiliary dual node id 2 lexicographically
    auxiliary_onprimaledge_dualedges::SVector{2, SVector{2, Any}}
    auxiliary_onprimalface_dualedges::SVector{2, SVector{2, Any}}
    Auxiliary_dualfacestruct() = new()
end


# declare struct to contain dual volume information
# a dual volume corresponds to a primal node
mutable struct Dualvolumestruct
    id::Int # primal node id

    interior_dualnodes::Vector{Int}
    boundary_dualnodes::Vector{SVector{3, Int}}
    auxiliary_dualnodes::Vector{SVector{2, Int}}

    interior_dualedges::Vector{SVector{3, Int}}
    boundary_dualedges::Vector{SVector{3, Int}}
    auxiliary_onprimaledge_dualedges::Vector{SVector{2, Any}}
    auxiliary_onprimalface_dualedges::Vector{SVector{2, Any}}

    interior_dualfaces::Vector{SVector{2, Int}}
    boundary_dualfaces::Vector{SVector{2, Int}}
    auxiliary_dualfaces::Vector{SVector{2, Any}}
    
    raw_volume::Float64
    Dualvolumestruct() = new()
end


# declare struct to contain all dual node dicts
mutable struct Dualnodedicts_struct
    interior_dualnodedict::Dict{Int, Interior_dualnodestruct}
    boundary_dualnodedict::Dict{SVector{3, Int}, Boundary_dualnodestruct}
    auxiliary_dualnodedict::Dict{SVector{2, Int}, Auxiliary_dualnodestruct}
    Dualnodedicts_struct() = new()
end


# declare struct to contain all dual edge dicts
mutable struct Dualedgedicts_struct
    interior_dualedgedict::Dict{SVector{3, Int}, Interior_dualedgestruct}
    boundary_dualedgedict::Dict{SVector{3, Int}, Boundary_dualedgestruct}
    auxiliary_onprimaledge_dualedgedict::Dict{SVector{2, Any}, Auxiliary_onprimaledge_dualedgestruct}
    auxiliary_onprimalface_dualedgedict::Dict{SVector{2, Any}, Auxiliary_onprimalface_dualedgestruct}
    Dualedgedicts_struct() = new()
end


# declare struct to contain all dual face dicts
mutable struct Dualfacedicts_struct
    interior_dualfacedict::Dict{SVector{2, Int}, Interior_dualfacestruct}
    boundary_dualfacedict::Dict{SVector{2, Int}, Boundary_dualfacestruct}
    auxiliary_dualfacedict::Dict{SVector{2, Any}, Auxiliary_dualfacestruct}
    Dualfacedicts_struct() = new()
end


# declare struct to contain all dual mesh dicts
mutable struct Dualmeshstruct
    dualnodedicts::Dualnodedicts_struct
    dualedgedicts::Dualedgedicts_struct
    dualfacedicts::Dualfacedicts_struct
    dualvolumedict::Dict{Int, Dualvolumestruct}
    Dualmeshstruct() = new()
end


########################################### End dual mesh ###########################################