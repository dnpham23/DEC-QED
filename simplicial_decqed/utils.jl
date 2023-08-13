"""
    vboundary(xmax::Float64,xmin::Float64,ymax::Float64,ymin::Float64,zmax::Float64,zmin::Float64, primalmesh::Primalmeshstruct, tol::Float64)

find boundary nodes of a 3D rectangular boundary
"""
function vboundary(xmax::Float64,xmin::Float64,ymax::Float64,ymin::Float64,zmax::Float64,zmin::Float64, primalmesh::Primalmeshstruct, tol::Float64)
    
    vbound = []
    
    for k in keys(primalmesh.nodedict)
        x = primalmesh.nodedict[k].coords[1]
        y = primalmesh.nodedict[k].coords[2]
        z = primalmesh.nodedict[k].coords[3]
        
        if (abs(x-xmax)<=tol)||(abs(x-xmin)<=tol)||(abs(y-ymax)<=tol)||(abs(y-ymin)<=tol)||(abs(z-zmax)<=tol)||(abs(z-zmin)<=tol)
            push!(vbound,k)
        end
    end
    
    return vbound
end

"""
    eboundary(vbound::Array{Any,1}, primalmesh::Primalmeshstruct, tol::Float64)

find boundary edges
"""
function eboundary(vbound::Array{Any,1}, primalmesh::Primalmeshstruct, tol::Float64)
    ebound = []
    
    for k in keys(primalmesh.edgedict)
        v1 = k[1]
        v2 = k[2]
        x1 = primalmesh.nodedict[v1].coords[1]
        y1 = primalmesh.nodedict[v1].coords[2]
        z1 = primalmesh.nodedict[v1].coords[3]
        x2 = primalmesh.nodedict[v2].coords[1]
        y2 = primalmesh.nodedict[v2].coords[2]
        z2 = primalmesh.nodedict[v2].coords[3]
        if (v1 in vbound)&&(v2 in vbound)&&( (abs(x1-x2)<=tol)||(abs(y1-y2)<=tol)||(abs(z1-z2)<=tol) )
            push!(ebound, k)
        end
    end
    return ebound
end


"""
    eboundary_gm(primalmesh::Primalmeshstruct, physicalnames_dict::Dict{Int, Physicalname_struct}, all_entities_struct::All_entities_struct)

find boundary edges
"""
function eboundary_gm(edgedict::Dict{SVector{2, Int}, Edgestruct}, physicalnames_dict::Dict{Int, Physicalname_struct}, all_entities_struct::All_entities_struct)
    ebound = []
    physicalnames_keys = collect(keys(physicalnames_dict)); 
    boundary_2dphysicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary1\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary2\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary3\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary4\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary5\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary6\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary7\"",  physicalnames_keys)][1],
                               physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"boundary8\"",  physicalnames_keys)][1]];
    # find the 2d entityid's associated with the boundaries
    surface_entities_keys = collect(keys(all_entities_struct.surface_entities_dict));
    bound_2d_entityid = surface_entities_keys[findall(x -> findall(in(all_entities_struct.surface_entities_dict[x].physicaltags), boundary_2dphysicaltags)!=[], surface_entities_keys)];
    for ekey in keys(primalmesh.edgedict)
        for ent_id in bound_2d_entityid
            if ent_id in edgedict[ekey].entities_dict[2] # if the current edge belongs to a superconducting entity
                if !(ekey in ebound)
                    push!(ebound, ekey)
                end
            end
        end
    end
    return ebound
end


"""
    insideTetCheck(Point::SArray{Tuple{3},Float64,1,3}, Tet::Array{SArray{Tuple{3},Float64,1,3},1})

    check whether a point lies inside a tetrahedron
"""
function insideTetCheck(Point::SArray{Tuple{3},Float64,1,3}, Tet::Array{SArray{Tuple{3},Float64,1,3},1})
    V1 = Tet[1]
    V2 = Tet[2]
    V3 = Tet[3]
    V4 = Tet[4]

    D0 = det([V1[1] V1[2] V1[3] 1.0; V2[1] V2[2] V2[3] 1.0; V3[1] V3[2] V3[3] 1.0; V4[1] V4[2] V4[3] 1.0])
    D1 = det([Point[1] Point[2] Point[3] 1.0; V2[1] V2[2] V2[3] 1.0; V3[1] V3[2] V3[3] 1.0; V4[1] V4[2] V4[3] 1.0])
    D2 = det([V1[1] V1[2] V1[3] 1.0; Point[1] Point[2] Point[3] 1.0; V3[1] V3[2] V3[3] 1.0; V4[1] V4[2] V4[3] 1.0])
    D3 = det([V1[1] V1[2] V1[3] 1.0; V2[1] V2[2] V2[3] 1.0; Point[1] Point[2] Point[3] 1.0; V4[1] V4[2] V4[3] 1.0]) 
    D4 = det([V1[1] V1[2] V1[3] 1.0; V2[1] V2[2] V2[3] 1.0; V3[1] V3[2] V3[3] 1.0; Point[1] Point[2] Point[3] 1.0])

    if (sign(D0)!=sign(D1))||(sign(D0)!=sign(D2))||(sign(D0)!=sign(D3))||(sign(D0)!=sign(D4))
        relpos = 0
    else
        relpos = 1
    end
    return relpos
end


"""
    insideTriCheck(Point::SArray{Tuple{3},Float64,1,3}, Tri::Array{SArray{Tuple{3},Float64,1,3},1})

    check whether a point lies inside a triangle
"""
function insideTriCheck(Point::SArray{Tuple{3},Float64,1,3}, Tri::Array{SArray{Tuple{3},Float64,1,3},1})
    V1 = Tri[1]
    V2 = Tri[2]
    V3 = Tri[3]
   
    # choose V1 as the origin of the barycentric coords
    vec0 = V2-V1
    vec1 = V3-V1
    vecP = Point-V1

    D00 = dot(vec0,vec0)
    D01 = dot(vec0,vec1)
    D0p = dot(vec0,vecP)
    D11 = dot(vec1,vec1)
    D1p = dot(vec1,vecP)

    # find the barycentric coordinates of Point
    DD = D00*D11 - D01^2
    u  = (D11*D0p - D01*D1p)/DD
    v  = (D00*D1p - D01*D0p)/DD

    if ((u>=0)&&(v>=0)&&(u+v<1))
        relpos = 1
    else 
        relpos = 0
    end
    return relpos
end


"""
    effective_material(nodedict::Dict{Int, Nodestruct}, 
                            edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                            facedict::Dict{SVector{3, Int}, Facestruct}, 
                            tetdict::Dict{Int, Tetstruct}, 
                            material::Array{Float64,1}, 
                            dualfacedicts::Dualfacedicts_struct)::Array{Float64,1}

compute the effective value on primal edges of a material property
"""
function effective_material(nodedict::Dict{Int, Nodestruct}, 
                            edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                            facedict::Dict{SVector{3, Int}, Facestruct}, 
                            tetdict::Dict{Int, Tetstruct}, 
                            material::Array{Float64,1}, 
                            dualfacedicts::Dualfacedicts_struct)::Array{Float64,1}

    interioredge_key = keys(dualfacedicts.interior_dualfacedict)
    boundaryedge_key = keys(dualfacedicts.boundary_dualfacedict) 

    material_eff = []
    for ekey in keys(edgedict)
        if ekey in interioredge_key
            append!(material_eff, get_dualarea_effectivevalue(edgedict[ekey], nodedict, facedict, tetdict, material)/dualfacedicts.interior_dualfacedict[ekey].raw_area)
        elseif ekey in boundaryedge_key
            append!(material_eff, get_dualarea_effectivevalue(edgedict[ekey], nodedict, facedict, tetdict, material)/dualfacedicts.boundary_dualfacedict[ekey].raw_area)
        end
    end

    return material_eff
end


"""
    get_effective_supportvolume(nodedict::Dict{Int, Nodestruct}, 
                                    edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                                    facedict::Dict{SVector{3, Int}, Facestruct}, 
                                    tetdict::Dict{Int, Tetstruct}, 
                                    material::Array{Float64,1}, 
                                    dualfacedicts::Dualfacedicts_struct)::Array{Float64,1}

compute the effective support volumes of primal edges
"""
function get_effective_supportvolume(nodedict::Dict{Int, Nodestruct}, 
                                    edgedict::Dict{SVector{2, Int}, Edgestruct}, 
                                    facedict::Dict{SVector{3, Int}, Facestruct}, 
                                    tetdict::Dict{Int, Tetstruct}, 
                                    material::Array{Float64,1}, 
                                    dualfacedicts::Dualfacedicts_struct)::Array{Float64,1}

    interioredge_key = keys(dualfacedicts.interior_dualfacedict)
    boundaryedge_key = keys(dualfacedicts.boundary_dualfacedict) 

    supportvolume_eff = [] 
    for ekey in keys(edgedict)
        if ekey in interioredge_key
            append!(supportvolume_eff, get_dualarea_effectivevalue(edgedict[ekey], nodedict, facedict, tetdict, material)*edgedict[ekey].length)
        elseif ekey in boundaryedge_key
            append!(supportvolume_eff, get_dualarea_effectivevalue(edgedict[ekey], nodedict, facedict, tetdict, material)*edgedict[ekey].length)
        end
    end

    return supportvolume_eff
end


"""
    triAngleCompute()

    find the angles of a triangle
"""
function triAngleCompute(Point1::Array{Float64,1}, Point2::Array{Float64,1}, Point3::Array{Float64,1})
    L12 = norm(Point1-Point2)
    L23 = norm(Point2-Point3)
    L31 = norm(Point3-Point1)

    A1 = acos((L12^2+L31^2-L23^2)/(2*L12*L31))
    A2 = acos((L12^2+L23^2-L31^2)/(2*L12*L23))
    A3 = acos((L23^2+L31^2-L12^2)/(2*L23*L31))

    return A1, A2, A3
end


"""
    halfspace_line_check()

    check whether two points lie in the same hal-space created by a line
"""
###function halfspace_line_check(Point1::Array{Float64,1}, Point2::Array{Float64,1}, Line::Array{SArray{Tuple{3},Float64,1,3},1})
function halfspace_line_check(Point1::SArray{Tuple{3},Float64,1,3}, Point2::SArray{Tuple{3},Float64,1,3}, Line::Array{SArray{Tuple{3},Float64,1,3},1})

    Line    = Array{Array{Float64,1},1}(Line)
    linevec = Line[2] - Line[1]
    Point1  = Array{Float64, 1}(Point1)
    Point2  = Array{Float64, 1}(Point2)
    ###return isequal(sign(cross(Point1-Line[1], linevec)),sign(cross(Point2-Line[1], linevec)))
    return isequal(sign(dot(cross(Point1-Line[1], linevec), cross(Point2-Line[1], linevec))), 1)
end


"""
    findSphericalCoords(x::Float64, y::Float64, z::Float64)

    find the spherical coordinates of a point given its cartesian coordinates
"""
function findSphericalCoords(x::Float64, y::Float64, z::Float64)
    r     = sqrt(x^2 + y^2 + z^2)
    theta = acos(z/r)
    if (x>0)&&(y>=0) # 1st quadrant
        omega = atan(y/x)
    elseif (x<=0)&&(y>0) # second quadrant
        omega = pi-atan(-y/x)
    elseif (x<0)&&(y<=0) # third quadrant
        omega = pi+atan(y/x)
    elseif (x>=0)&&(y<0) # 4th quadrant
        omega = 2*pi-atan(-y/x)
    end
    return r, theta, omega
end


"""
    find_wg_edges(edgedict::Dict{SVector{2, Int}, Edgestruct},
                       nodedict::Dict{Int, Nodestruct},
                       physicalnames_dict::Dict{Int, Physicalname_struct}, 
                       all_entities_struct::All_entities_struct, 
                       t::Float64, tol::Float64)

find the edges that belong to superconducting entitites
"""
function find_sc_edges(edgedict::Dict{SVector{2, Int}, Edgestruct},
                       nodedict::Dict{Int, Nodestruct},
                       physicalnames_dict::Dict{Int, Physicalname_struct}, 
                       all_entities_struct::All_entities_struct)

    sc_edgelist  = []
    # find the physical tags associated with superconducting objects
    physicalnames_keys = collect(keys(physicalnames_dict)); 
    sc_2dphysicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor1\"",  physicalnames_keys)][1],
                         physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor2\"",  physicalnames_keys)][1],
                         physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor3\"",  physicalnames_keys)][1],
                         physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor4\"",  physicalnames_keys)][1]];
    sc_1dphysicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"JJ1\"",  physicalnames_keys)][1], 
                         physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"JJ2\"",  physicalnames_keys)][1]];
    # find the 2d entityid's associated with superconducting objects
    surface_entities_keys = collect(keys(all_entities_struct.surface_entities_dict));
    sc_2d_entityid = surface_entities_keys[findall(x -> findall(in(all_entities_struct.surface_entities_dict[x].physicaltags), sc_2dphysicaltags)!=[], surface_entities_keys)];
    # find the 1d entityid's associated with superconducting objects
    curve_entities_keys = collect(keys(all_entities_struct.curve_entities_dict));
    sc_1d_entityid = curve_entities_keys[findall(x -> findall(in(all_entities_struct.curve_entities_dict[x].physicaltags), sc_1dphysicaltags)!=[], curve_entities_keys)];
    for ekey in keys(edgedict)
        for ent_id in sc_2d_entityid
            if ent_id in edgedict[ekey].entities_dict[2] # if the current edge belongs to a superconducting entity
                if !(ekey in sc_edgelist)
                    push!(sc_edgelist, ekey)
                end
            end
        end

        for ent_id in sc_1d_entityid
            if ent_id in edgedict[ekey].entities_dict[1] # if the current edge belongs to a superconducting entity
                if !(ekey in sc_edgelist)
                    push!(sc_edgelist, ekey)
                end
            end 
        end
    end

    return sc_edgelist;
end


"""
    find_substr_edges(edgedict::Dict{SVector{2, Int}, Edgestruct},
                           physicalnames_dict::Dict{Int, Physicalname_struct}, 
                           all_entities_struct::All_entities_struct, 
                           sc_edgelist::Array{Any,1})

find edges that belong to the substrates
"""
function find_substr_edges(edgedict::Dict{SVector{2, Int}, Edgestruct},
                           physicalnames_dict::Dict{Int, Physicalname_struct}, 
                           all_entities_struct::All_entities_struct, 
                           sc_edgelist::Array{Any,1})
    substr_edgelist = []
    # find the physical tags associated with the substrate
    physicalnames_keys  = collect(keys(physicalnames_dict)); 
    substr_physicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"substrate1\"",  physicalnames_keys)][1],
                           physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"substrate2\"",  physicalnames_keys)][1]];
    # find the 3d entityid's associated with the substrate
    volume_entities_keys = collect(keys(all_entities_struct.volume_entities_dict));
    substr_3d_entityid = volume_entities_keys[findall(x -> findall(in(all_entities_struct.volume_entities_dict[x].physicaltags), substr_physicaltags)!=[], volume_entities_keys)];
    for ekey in keys(edgedict)
        for ent_id in substr_3d_entityid
            if ent_id in edgedict[ekey].entities_dict[3]
                # need to eliminate edges lying on the surface of the substrate that belongs to the superconducting components
                if !(ekey in sc_edgelist)
                    if !(ekey in substr_edgelist)
                        push!(substr_edgelist, ekey)
                    end
                end
            end
        end
    end

    return substr_edgelist;
end


function find_substr1_surface_edge(edgedict::Dict{SVector{2, Int}, Edgestruct},
                           physicalnames_dict::Dict{Int, Physicalname_struct}, 
                           all_entities_struct::All_entities_struct)

    substr1_surface_edgelist = []
    # find the physical tags associated with the substrate
    physicalnames_keys  = collect(keys(physicalnames_dict)); 
    substr1_surface_physicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"substrate 1's surface\"",  physicalnames_keys)][1],
                                    physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor1\"",  physicalnames_keys)][1],
                                    physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"superconductor2\"",  physicalnames_keys)][1]];
    JJ1_physicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"JJ1\"",  physicalnames_keys)][1]];
    # find the 2d entityid's associated with the substrate
    surface_entities_keys = collect(keys(all_entities_struct.surface_entities_dict));
    substr1_2d_entityid = surface_entities_keys[findall(x -> findall(in(all_entities_struct.surface_entities_dict[x].physicaltags), substr1_surface_physicaltags)!=[], surface_entities_keys)];
    # find the 1d entityid's associated with JJ1
    curve_entities_keys = collect(keys(all_entities_struct.curve_entities_dict));
    JJ1_1d_entityid = curve_entities_keys[findall(x -> findall(in(all_entities_struct.curve_entities_dict[x].physicaltags), JJ1_physicaltags)!=[], curve_entities_keys)];

    for ekey in keys(edgedict)
        for ent_id in substr1_2d_entityid
            if ent_id in edgedict[ekey].entities_dict[2] # if the current edge belongs to a superconducting entity
                if !(ekey in substr1_surface_edgelist)
                    push!(substr1_surface_edgelist, ekey)
                end
            end
        end

        for ent_id in JJ1_1d_entityid
            if ent_id in edgedict[ekey].entities_dict[1] # if the current edge belongs to a superconducting entity
                if !(ekey in substr1_surface_edgelist)
                    push!(substr1_surface_edgelist, ekey)
                end
            end 
        end
    end

    return substr1_surface_edgelist;
end

function find_JJ1_nodes(nodedict::Dict{Int, Nodestruct},
                        physicalnames_dict::Dict{Int, Physicalname_struct}, 
                        all_entities_struct::All_entities_struct)
    JJ1_nodelist = []
    # find the physical tags associated with the substrate
    physicalnames_keys  = collect(keys(physicalnames_dict)); 
    JJ1_physicaltags = [physicalnames_keys[findall(x -> physicalnames_dict[x].name == "\"JJ1\"",  physicalnames_keys)][1]];
    # find the 1d entityid's associated with JJ1
    curve_entities_keys = collect(keys(all_entities_struct.curve_entities_dict));
    JJ1_1d_entityid = curve_entities_keys[findall(x -> findall(in(all_entities_struct.curve_entities_dict[x].physicaltags), JJ1_physicaltags)!=[], curve_entities_keys)];

    for nkey in keys(nodedict)
        for ent_id in JJ1_1d_entityid
            if ent_id in nodedict[nkey].entities_dict[1] # if the current edge belongs to a superconducting entity
                if !(nkey in JJ1_nodelist)
                    push!(JJ1_nodelist, nkey)
                end
            end 
        end
    end
    return JJ1_nodelist;
end