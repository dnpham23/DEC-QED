# This file is used to implement the functionality to parse the .msh file.
# Note that information about the material imhomogeneities are captured
# by the parsephysicalnames and parse_entities functions, but they do 
# not currently used in `primalmesh.jl` or `dualmesh.jl``.

include("types.jl")


"""
    getfilepath(path::String)::String

Return the string of the complete path to file.

This function concatenates the path referencing the parent directory to the path of file.
"""
function getfilepath(file::String)::String

    parentdir = string(pwd())
    filepath = parentdir * file

end


"""
    parsephysicalnames(fileinstance::IOStream)

Reads fileinstance line by line.
"""
function parsephysicalnames(fileinstance::IOStream)::Dict{Int, Physicalname_struct}

    # instantiate empty dict to populate 
    physicalnames_dict = Dict{Int, Physicalname_struct}()

    # 1st line after $PhysicalNames gives number of physicalnames
    currentline = readline(fileinstance)
    num_physicalnames = parse.(Int, split(strip(currentline), " "))[1]
    
    for _ in 1:num_physicalnames

        currentline = readline(fileinstance) 
        physicalname_data = split(strip(currentline), " ")

        dimension = parse(Int, physicalname_data[1])
        physicaltag = parse(Int, physicalname_data[2])
        name = join(physicalname_data[3:end], " ")

        physicalname_struct = Physicalname_struct()
        physicalname_struct.physicaltag = physicaltag
        physicalname_struct.dimension = dimension
        physicalname_struct.name = name

        physicalnames_dict[physicaltag] = physicalname_struct

    end

    return physicalnames_dict

end


"""

Return an dict of entities corresponding to the relevant
information of the "\$Entities" section in fileinstance.
"""
function parse_entities_2(fileinstance::IOStream)::All_entities_struct

    # instantiate empty struct to populate
    all_entities_struct = All_entities_struct()

    # 1st line after $Entities gives number of entities of each category
    currentline = readline(fileinstance)
    entities_data = parse.(Int, split(strip(currentline), " "))
    num_points = entities_data[1]
    num_curves = entities_data[2]
    num_surfaces = entities_data[3]
    num_volumes = entities_data[4]

    # create dictionary of point entities
    point_entities_dict = Dict{Int, Point_entity_struct}()
    for _ in 1:num_points

        currentline = readline(fileinstance)
        entitydata = split(strip(currentline), " ")
        point_tag = parse(Int, entitydata[1])
        num_physicaltags = parse(Int, entitydata[5])
        
        physicaltags = []
        if num_physicaltags != 0
            append!(physicaltags, parse.(Int, entitydata[6:6+num_physicaltags-1]))
        end

        point_entity_struct = Point_entity_struct()
        point_entity_struct.point_tag = point_tag
        point_entity_struct.physicaltags = physicaltags

        point_entities_dict[point_tag] = point_entity_struct

    end

    # create dictionary of curve entities
    curve_entities_dict = Dict{Int, Curve_entity_struct}()
    for _ in 1:num_curves

        currentline = readline(fileinstance)
        entitydata = split(strip(currentline), " ")
        curve_tag = parse(Int, entitydata[1])
        num_physicaltags = parse(Int, entitydata[8])
        num_boundingPoints = parse(Int, entitydata[9+num_physicaltags])

        physicaltags = []
        if num_physicaltags != 0
            append!(physicaltags, parse.(Int, entitydata[9:9+num_physicaltags-1]))
        end

        boundingPoints = []
        bnd_pts_startind = 9+num_physicaltags+1
        append!(boundingPoints, abs.(parse.(Int, entitydata[bnd_pts_startind:bnd_pts_startind+num_boundingPoints-1])))

        curve_entity_struct = Curve_entity_struct()
        curve_entity_struct.curve_tag = curve_tag
        curve_entity_struct.physicaltags = physicaltags
        curve_entity_struct.boundingpoints = boundingPoints

        curve_entities_dict[curve_tag] = curve_entity_struct

    end

    # create dictionary of surface entities
    surface_entities_dict = Dict{Int, Surface_entity_struct}()
    for _ in 1:num_surfaces

        currentline = readline(fileinstance)
        entitydata = split(strip(currentline), " ")
        surface_tag = parse(Int, entitydata[1])
        num_physicaltags = parse(Int, entitydata[8])
        num_boundingCurves = parse(Int, entitydata[9+num_physicaltags])

        physicaltags = []
        if num_physicaltags != 0
            append!(physicaltags, parse.(Int, entitydata[9:9+num_physicaltags-1]))
        end

        boundingCurves = []
        bnd_curves_startind = 9+num_physicaltags+1
        append!(boundingCurves, abs.(parse.(Int, entitydata[bnd_curves_startind:bnd_curves_startind+num_boundingCurves-1])))

        surface_entity_struct = Surface_entity_struct()
        surface_entity_struct.surface_tag = surface_tag
        surface_entity_struct.physicaltags = physicaltags
        surface_entity_struct.boundingcurves = boundingCurves

        surface_entities_dict[surface_tag] = surface_entity_struct

    end

    # create dictionary of volume entities
    volume_entities_dict = Dict{Int, Volume_entity_struct}()
    for _ in 1:num_volumes

        currentline = readline(fileinstance)
        entitydata = split(strip(currentline), " ")
        volume_tag = parse(Int, entitydata[1])
        num_physicaltags = parse(Int, entitydata[8])
        num_boundingSurfaces = parse(Int, entitydata[9+num_physicaltags])

        println(num_physicaltags)
        
        physicaltags = []
        if num_physicaltags != 0
            append!(physicaltags, parse.(Int, entitydata[9:9+num_physicaltags-1]))
        end

        boundingSurfaces = []
        bnd_surfaces_startind = 9+num_physicaltags+1
        append!(boundingSurfaces, abs.(parse.(Int, entitydata[bnd_surfaces_startind:bnd_surfaces_startind+num_boundingSurfaces-1])))

        volume_entity_struct = Volume_entity_struct()
        volume_entity_struct.volume_tag = volume_tag
        volume_entity_struct.physicaltags = physicaltags
        volume_entity_struct.boundingsurfaces = boundingSurfaces

        volume_entities_dict[volume_tag] = volume_entity_struct
    end

    # insert dicts into struct
    all_entities_struct.point_entities_dict = point_entities_dict
    all_entities_struct.curve_entities_dict = curve_entities_dict
    all_entities_struct.surface_entities_dict = surface_entities_dict
    all_entities_struct.volume_entities_dict = volume_entities_dict

    return all_entities_struct

end





"""
    parsenodes(fileinstance::IOStream)

Return a dict of Nodestruct corresponding to the "\$Nodes" section in fileinstance.

Reads fileinstance line by line.
"""
function parsenodes_2(fileinstance::IOStream, all_entities_struct::All_entities_struct)::Dict{Int64, Nodestruct}

    # instantiate empty dict to populate with structs
    nodedict = Dict{Int, Nodestruct}()
    
    # 1st line after $Nodes contains metadata about section
    currentline = readline(fileinstance)
    nodesdata = parse.(Int, split(strip(currentline), " "))
    num_node_entities = nodesdata[1]
    # num_nodes = nodesdata[2]

    for _ in 1:num_node_entities
        
        # 1st line of each block contains meta-data about entity
        currentline = readline(fileinstance) 
        entitydata = parse.(Int, split(strip(currentline), " "))
        root_entitydim = entitydata[1] # MshFileVersion 4.1
        root_entityid = entitydata[2] # index of the entity (MshFileVersion 4.1)
        num_nodes_inentity = entitydata[4]

        # list of all node tags in this entity
        current_tags = Vector{Int}() 
        
        # if entity has x nodes, first x lines after 
        # entity meta-data contains all node tags in this entity
        for _ in 1:num_nodes_inentity
            
            currentline = readline(fileinstance)
            tag = parse.(Int, currentline)
            push!(current_tags, tag)
        
        end

        for nodeid in current_tags

            currentline = readline(fileinstance)
            nodecoords = parse.(Float64, split(strip(currentline), " ")) 
            nodecoords = convert(SVector{3, Float64}, nodecoords)

            entities_dict = Dict{Int,Vector{}}()
            # from the root entity, find all higher-dim entities that this node belongs to
            if root_entitydim == 0
                curvetags = []
                for kcurve in keys(all_entities_struct.curve_entities_dict)
                    if root_entityid in all_entities_struct.curve_entities_dict[kcurve].boundingpoints
                        append!(curvetags, kcurve)
                    end
                end
                surfacetags = []
                for ksurface in keys(all_entities_struct.surface_entities_dict)
                    for kcurve in curvetags
                        if (kcurve in all_entities_struct.surface_entities_dict[ksurface].boundingcurves)&&!(ksurface in surfacetags)
                            append!(surfacetags, ksurface)
                        end
                    end
                end
                volumetags = []
                for kvolume in keys(all_entities_struct.volume_entities_dict)
                    for ksurface in surfacetags
                        if (ksurface in all_entities_struct.volume_entities_dict[kvolume].boundingsurfaces)&&!(kvolume in volumetags)
                            append!(volumetags, kvolume)
                        end
                    end
                end
                entities_dict[1] = curvetags
                entities_dict[2] = surfacetags
                entities_dict[3] = volumetags
            elseif root_entitydim == 1
                surfacetags = []
                for ksurface in keys(all_entities_struct.surface_entities_dict)
                    if (root_entityid in all_entities_struct.surface_entities_dict[ksurface].boundingcurves)&&!(ksurface in surfacetags)
                        append!(surfacetags, ksurface)
                    end
                end
                volumetags = []
                for kvolume in keys(all_entities_struct.volume_entities_dict)
                    for ksurface in surfacetags
                        if (ksurface in all_entities_struct.volume_entities_dict[kvolume].boundingsurfaces)&&!(kvolume in volumetags)
                            append!(volumetags, kvolume)
                        end
                    end
                end
                entities_dict[1] = [root_entityid]
                entities_dict[2] = surfacetags
                entities_dict[3] = volumetags
            elseif root_entitydim == 2
                volumetags = []
                for kvolume in keys(all_entities_struct.volume_entities_dict)
                    if (root_entityid in all_entities_struct.volume_entities_dict[kvolume].boundingsurfaces)&&!(kvolume in volumetags)
                        append!(volumetags, kvolume)
                    end
                end
                entities_dict[1] = [0]
                entities_dict[2] = [root_entityid]
                entities_dict[3] = volumetags
            elseif root_entitydim == 3
                entities_dict[1] = [0]
                entities_dict[2] = [0]
                entities_dict[3] = [root_entityid]
            end

            # create a Nodestruct
            node = Nodestruct()
            node.id = nodeid
            node.coords = nodecoords
            node.root_entitydim = root_entitydim
            node.root_entityid  = root_entityid
            node.entities_dict  = entities_dict
            # update nodedict
            nodedict[nodeid] = node

        end 

    end

    return nodedict

end


"""
    parsetets(nodedict::Vector{Nodestruct}, fileinstance::IOStream)

Return a dict of structs of Tetstruct corresponding to the "\$Elements" section in fileinstance.

Reads fileinstance line by line.
The output is a Vector of Tetstruct, which are each mutable. 
"""
function parsetets(fileinstance::IOStream)::Dict{Int, Tetstruct}

    # instantiate empty dict to populate with structs
    tetdict = Dict{Int, Tetstruct}()

    # 1st line after $elements contains metadata about section
    currentline = readline(fileinstance)
    elementsdata = parse.(Int, split(strip(currentline), " "))
    num_element_entities = elementsdata[1]
    num_elements_inentity = elementsdata[2]

    for _ in 1:num_element_entities
        
        # 1st line of each block contains meta-data about entity
        currentline = readline(fileinstance) 
        entitydata = parse.(Int, split(strip(currentline), " "))
        entityid = entitydata[2]
        elementtype = entitydata[3]
        num_elements_inentity = entitydata[4]

        if elementtype == 4
            
            for _ in 1:num_elements_inentity

                currentline = readline(fileinstance)
                currentline_parsed = parse.(Int, split(strip(currentline), " "))
                elementtag = currentline_parsed[1]
                element_nodeids = currentline_parsed[2:end] # 4 Ints

                # create a Tetstruct
                tet = Tetstruct()

                tet.id = elementtag
                tet.nodes = element_nodeids
                tet.entityid = entityid

                # update tetdict
                tetdict[elementtag] = tet
                
            end

        else
            for _ in 1:num_elements_inentity
                currentline = readline(fileinstance)
            end
        end
    
    end

    return tetdict

end 


"""
    parsefile(file::String)

Return the dictionaries of nodes, tets, physical names, and entities.

Reads the .msh file line by line and calls the appropriate section
parser.
"""
function parsefile(file::String)::Vector{Any}

    physicalnames_dict = Dict{Int, Physicalname_struct}()
    all_entities_struct = All_entities_struct()
    nodedict = Dict{Int, Nodestruct}()
    tetdict = Dict{Int, Nodestruct}()

    # get full file path of .msh file
    meshpath = getfilepath(file)

    # parse file
    currentline = "init" 
    open(meshpath) do fileinstance

        while !eof(fileinstance)

            if startswith(currentline, "\$PhysicalNames")
                physicalnames_dict = parsephysicalnames(fileinstance)
            elseif startswith(currentline, "\$Entities")
                all_entities_struct = parse_entities_2(fileinstance)
            elseif startswith(currentline, "\$Nodes")
                nodedict = parsenodes_2(fileinstance, all_entities_struct)
            elseif startswith(currentline, "\$Elements")
                tetdict = parsetets(fileinstance)
            end
            currentline = readline(fileinstance)

        end

    end

    return [nodedict, tetdict, physicalnames_dict, all_entities_struct]

end



