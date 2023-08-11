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