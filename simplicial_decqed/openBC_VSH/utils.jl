"""
    doublecurl_openbc(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound3::Array{Any,1}, ebound_extra::Array{SArray{Tuple{2},Int64,1,2},1},
                           edgedict_extra::Dict{SArray{Tuple{2},Int64,1,2},Edgestruct}, edgekeylist::Array{SArray{Tuple{2},Int64,1,2},1})

compute curl-curl operator for VHS boundary condition
"""
function doublecurl_openbc(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound3::Array{Any,1}, ebound_extra::Array{SArray{Tuple{2},Int64,1,2},1},
                           edgedict_extra::Dict{SArray{Tuple{2},Int64,1,2},Edgestruct}, edgekeylist::Array{SArray{Tuple{2},Int64,1,2},1})
    Ne = length(primalmesh.edgedict) + length(edgedict_extra)
    dcurl = Array{Float64,2}(undef, Ne, Ne)
    
    facekeylist = collect(keys(primalmesh.facedict))

    for i = 1:Ne
        ki = edgekeylist[i]
        if !(ki in ebound3)&&!(ki in ebound_extra)
            v1 = ki[1]
            v2 = ki[2]
            dual_area  = dualmesh.dualfacedicts.interior_dualfacedict[ki].raw_area
            primal_len = primalmesh.edgedict[ki].length
            for j = 1:Ne
                kj = edgekeylist[j]
                if (v1==kj[1])&&(v2!=kj[2]) # if the current edge and neighboring edge share v1 node
                    facekey       = sort([v1,v2,kj[2]]) 
                    if (facekey in facekeylist)
                        primal_len_kj = primalmesh.edgedict[kj].length
                        area        = primalmesh.facedict[facekey].area
                        dual_len    = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j]  = -(dual_len/area)*(primal_len_kj/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop       
                        dcurl[i,i] += (dual_len/area)*(primal_len/dual_area)  # for each loop coming out of v1, add up contribution of the central edge
                    end
                elseif (v1==kj[2])&&(v2!=kj[1]) # if the current edge and neighboring edge share v2 node
                    facekey       = sort([v1,v2,kj[1]]) 
                    if (facekey in facekeylist)
                        primal_len_kj = primalmesh.edgedict[kj].length
                        area        = primalmesh.facedict[facekey].area
                        dual_len    = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j]  = (dual_len/area)*(primal_len_kj/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                        dcurl[i,i] += (dual_len/area)*(primal_len/dual_area) # for each loop coming out of v1, add up contribution of the central edge
                    end
                elseif (v2==kj[1])&&(v1!=kj[2])
                    facekey       = sort([v1,v2,kj[2]]) 
                    if (facekey in facekeylist)
                        primal_len_kj = primalmesh.edgedict[kj].length
                        area       = primalmesh.facedict[facekey].area
                        dual_len   = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j] = (dual_len/area)*(primal_len_kj/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                    end
                elseif (v2==kj[2])&&(v1!=kj[1])
                    facekey       = sort([v1,v2,kj[1]]) 
                    if (facekey in facekeylist)
                        primal_len_kj = primalmesh.edgedict[kj].length
                        area       = primalmesh.facedict[facekey].area
                        dual_len   = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j] = -(dual_len/area)*(primal_len_kj/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                    end
                end
            end
        end
    end
    
    return dcurl
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
###function insideTriCheck(Point::Array{Float64,1}, Tri::Array{SArray{Tuple{3},Float64,1,3},1})
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
        relpos = true
    else 
        relpos = false
    end
    return relpos
end


# compute the primal-dual dot product between a primal 1-form and a dual 2-form. The result of this product is dual 3-form, which is converted to primal 0-form living on primal vertices
### NOT TESTED YET ###
function pd_dot(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, p1f::Array{Float64,2}, d2f::Array{Float64,2})
    Ne   = length(primalmesh.edgedict)
    Nv   = length(primalmesh.nodedict)
    nodekeylist = collect(keys(primalmesh.nodedict))
    edgekeylist = collect(keys(primalmesh.edgedict))
    dp_prod = []

    for iv = 1:Nv
        nodeid  = nodekeylist[iv]
        # find all primal edges attached to this vertex
        edgelist = []
        for edge in primalmesh.edgedict
            edgeid = edge.first
            if nodeid in edgeid
                push!(edgelist, edgeid)
            end
        end
        
        temp_prod = 0.0
        for edgeid in edgelist
            # find location of the current edge in the list of 1-form fields
            edge_ind = findfirst(x->x==edgeid, edgekeylist)
            temp_prod += 0.5*p1f[edge_ind]*d2f[edge_ind]  # this 0.5 factor is likely false (due to a typo in Hirani's definition of dot product defined at vertices). 
                                                          # should be the ratio (support volume ∩ dual volume)/(dual volume) instead
        end
	dual_vol = dualmesh.dualvolumedict[nodeid].raw_volume
	push!(dp_prod, temp_prod./dual_vol)
    end

    return dp_prod
end

# compute the quantum pressure term (1-form)
### NOT TESTED YET ###
function q_pressure(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, rho::Array{Float64,2})
    Ne   = length(primalmesh.edgedict)
    nodekeylist = collect(keys(primalmesh.nodedict))
    edgekeylist = collect(keys(primalmesh.edgedict))
    qp = []
    
    psi_amp = sqrt.(abs.(real(rho)))
    for ie = 1:Ne
        ki = edgekeylist[ie]
	if (ki in sc_edgelist) # if the edge is in a superconducting region
	    temp_lap = zeros(2)
            for iv = 1:2
	       v = ki[iv]
	       # find all primal edges attached to this vertex
               edgelist = []
               for edge in primalmesh.edgedict
                   edgeid = edge.first
                   if nodeid in edgeid
                       push!(edgelist, edgeid)
                   end
	       end

	       # find the location of the node in the list of 0-forms
               v_ind = findfirst(x->x==v, nodekeylist)
               for edgeid in edgelist
                   dual_area  = dualmesh.dualfacedicts.interior_dualfacedict[edgeid].raw_area
                   primal_len = primalmesh.edgedict[edgeid].length
	           other_v    = edgeid[findfirst(x->x!=v,edgeid)]
	           # find the location of the node in the list of 0-forms
	           other_v_ind   = findfirst(x->x==other_v, nodekeylist)
	           # accumulate value for the laplacian
	           temp_lap[iv] += dual_area/primal_len*(psi_amp[other_v_ind]-psi_amp[v_ind]) 
	       end
	       temp_lap[iv] = temp_lap[iv]./(psi_amp[v_ind]*dualmesh.dualvolumedict[v].raw_volume)
	    end
	    push!(qp,sign(ki[1]-ki[2])*(temp_lap[1]-temp_lap[2]))
	else
	    push!(qp,0.0)
        end
    end
    return qp # result is a set of 1-forms
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
    boundary_op(primalmesh::Primalmeshstruct, fbound2::Array{Any,1}, ebound2::Array{Any,1}, ebound3::Array{Any,1}, 
                    ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, ebound_extra::Array{Any,1}, vbound2::Array{Any,1},
                    r1::Float64, r2::Float64, r3::Float64, Lmax::Int64, k::Array{Any,1},
                    edgekeylist::Array{SArray{Tuple{2},Int64,1,2},1}, Ne::Int64)

compute the open-boundary operator
"""
function boundary_op(primalmesh::Primalmeshstruct, fbound2::Array{Any,1}, ebound2::Array{Any,1}, ebound3::Array{Any,1}, 
                    ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, ebound_extra::Array{SArray{Tuple{2},Int64,1,2},1}, vbound2::Array{Any,1},
                    r1::Float64, r2::Float64, r3::Float64, Lmax::Int64, Lmin::Int64, k::Complex{Float64},
                    edgekeylist::Array{SArray{Tuple{2},Int64,1,2},1}, Ne::Int64)

    bop = Array{Complex{Float64},2}(undef, Ne, Ne)

    for eb in ebound_extra
        eb_ind = findfirst(x->x==eb, edgekeylist)

        # compute the contributions of radial edges at inner boundary layers
        e_radial_coeff = Ar_coeff(primalmesh, eb, ebound_r12, ebound_r23, fbound2, vbound2, Lmax, Lmin, r1, r2, r3, k)
        e_radial_keys  = keys(e_radial_coeff)

        for ekey in e_radial_keys
            # find the position of this edge in the list of all edges
            ekey_ind = findfirst(x->x==ekey, edgekeylist)
            bop[eb_ind, ekey_ind] = e_radial_coeff[ekey]
        end
        bop[eb_ind, eb_ind] = -1.0
    end

    for eb in ebound3
        eb_ind = findfirst(x->x==eb, edgekeylist)

        # compute the contributions of tangential edges at inner boundary layer r2
        e_angle_coeff  = a_angle_coeff(primalmesh, eb, fbound2, ebound2, r2, r3, Lmax, Lmin, k)
        e_angle_keys  = keys(e_angle_coeff)

        for ekey in e_angle_keys
            # find the position of this edge in the list of all edges
            ekey_ind = findfirst(x->x==ekey, edgekeylist)
            bop[eb_ind, ekey_ind] = e_angle_coeff[ekey]
        end
        bop[eb_ind, eb_ind] = -1.0
    end


    return bop
end
