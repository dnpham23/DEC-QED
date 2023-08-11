using SpecialFunctions
using SphericalHarmonics

"""
    Ylm_normalizecheck(primalmesh::Primalmeshstruct, fbound1::Array{Any,1}, l::Int64, m::Int64 )

    check whether the Spherical Harmonics are normalized
"""
function Ylm_normalizecheck(primalmesh::Primalmeshstruct, fbound1::Array{Any,1}, l::Int64, m::Int64 )
    NormTemp = 0.0

    for fkey in fbound1
        circ = get_circumcenter_face(primalmesh.facedict[fkey], primalmesh.nodedict)
        x     = circ[1]
        y     = circ[2]
        z     = circ[3]
        r, theta, omega = findSphericalCoords(x, y, z)
        NormTemp += computeYlm(theta, omega, lmax = 2)[(2,-2)]*conj(computeYlm(theta, omega, lmax = 2)[(l,m)])*primalmesh.facedict[fkey].area
        #NormTemp += primalmesh.facedict[fkey].area
    end
    return NormTemp
end

"""
    Ar_coeff_full(primalmesh::Primalmeshstruct, edgeform::Array{Float64,1}, ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, 
                       fbound2::Array{Any,1}, Lm::Int64, r1::Float64, r2::Float64)

    compute the coefficients a_lm for open boundary conditions of Ar
"""
function Ar_coeff_full(primalmesh::Primalmeshstruct, edgeform::Array{Float64,1}, ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, 
                       fbound2::Array{Any,1}, Lm::Int64, r1::Float64, r2::Float64)
    A_lm_coeff = Dict{SVector{2, Int}, Float64}

    edegkeylist = collect(keys(edgedict))

    for l = 2:Lm # 0th and 1st order diverge at r = 0
        for m = -l:l
            akey = [l,m]
            temp_coeff = 0.0
            for fkey in fbound2
                meanform = 0.0
                # find the radial edges connected to the vertices of this face
                for i = 1:3
                    e1 = ebound_r12[findfirst(x->fkey[i] in x, ebound_r12)]
                    e2 = ebound_r23[findfirst(x->fkey[i] in x, ebound_r23)]
                    e1_ind = findfirst(x->x==e1, edegkeylist)
                    e2_ind = findfirst(x->x=e2, edegkeylist)
                    meanform += edgeform[e1_ind] + edgeform[e2_ind]
                end
                meanform = meanform/6

                circ = get_circumcenter_face(primalmesh.facedict[fkey], primalmesh.nodedict)
                r, theta, omega = findSphericalCoords(circ[1], circ[2], circ[3])
                temp_coeff += primalmesh.facedict[fkey].area*meanform*computeYlm(theta, omega, lmax = l)[(l,m)]
            end
            A_lm_coeff[akey] = 1/(sqrt(r2)*(r2-r1)*l*(l+1))*temp_coeff
        end
    end

    return A_lm_coeff
end

"""
    Ar_coeff(primalmesh::Primalmeshstruct, eb_current::Array{Any,1}, ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, 
                    fbound2::Array{Any,1}, vbound2::Array{Any,1}, Lm::Int64, r1::Float64, r2::Float64, r3::Float64, k)

    calculate the coefficients to be put in the matrix for open BC of the radial field (normal to the boundary sphere)
"""
function Ar_coeff(primalmesh::Primalmeshstruct, eb_current::SArray{Tuple{2},Int64,1,2}, ebound_r12::Array{Any,1}, ebound_r23::Array{Any,1}, 
                    fbound2::Array{Any,1}, vbound2::Array{Any,1}, Lmax::Int64, Lmin::Int64, r1::Float64, r2::Float64, r3::Float64, k::Complex{Float64})

    # the dictionary to store the coefficients (the keys are the edges that contribute to the boundary condition)
    ###eb_coeff = Dict{SVector{2, Int}, Float64}()
    eb_coeff = Dict{Array{Any, 1}, Complex{Float64}}()
    # find the spherical angles of the radial edge currently considered
    n1  = eb_current[1]
    #n2  = eb_current[2]
    x1  = primalmesh.nodedict[n1].coords[1]
    y1  = primalmesh.nodedict[n1].coords[2]
    z1  = primalmesh.nodedict[n1].coords[3]
    _, theta, phi = findSphericalCoords(x1, y1, z1)

    ebound_r = [ebound_r12; ebound_r23]
    dr23 = r3-r2
    dr12 = r2-r1
    for eb_key in ebound_r
        # find which node of this edge lies on the r2 sphere
        if eb_key[1] in vbound2
            # find the faces on fbound2 that this radial edge touches
            facelist = fbound2[findall(x->eb_key[1] in x, fbound2)]
        elseif eb_key[2] in vbound2
            facelist = fbound2[findall(x->eb_key[2] in x, fbound2)]
        end

        tempcoef = 0.0
        for fkey in facelist
            area = primalmesh.facedict[fkey].area
            circ = get_circumcenter_face(primalmesh.facedict[fkey], primalmesh.nodedict)
            _, theta_f, phi_f = findSphericalCoords(circ[1], circ[2], circ[3])
            for l = Lmin:Lmax
                for m = -l:l
                    tempcoef += (1 + dr23*( k*0.5*(hankelh1(l-1/2,k*r2)-hankelh1(l+3/2,k*r2))/hankelh1(l+1/2,k*r2) - 3/(2*r2) ) )*
                                    conj(computeYlm(theta_f, phi_f, lmax = l)[(l,m)])*computeYlm(theta, phi, lmax = l)[(l,m)]
                end
                ###tempcoef = tempcoef/(l*(l+1))
            end
            tempcoef = tempcoef*area
        end

        if eb_key in ebound_r12 
            if eb_key[1] in vbound2 # take care of the orientation of the 1-form \w 
                eb_coeff[eb_key] = -tempcoef*(2.0*dr23)/(6.0*dr12*r2^2)
            else
                eb_coeff[eb_key] = tempcoef*(2.0*dr23)/(6.0*dr12*r2^2)
            end
        else
            if eb_key[1] in vbound2 # take care of the orientation of the 1-form \w 
                eb_coeff[eb_key] = tempcoef*(2.0*dr23)/(6.0*dr23*r2^2)
            else
                eb_coeff[eb_key] = -tempcoef*(2.0*dr23)/(6.0*dr23*r2^2)
            end
        end
    end

    # contribution of the radial edge right beneath the current edge
    e_below_key = ebound_r23[findfirst(x->n1 in x, ebound_r23)]
    if e_below_key[1] == n1 # if this beneath edge points into the center (-r direction)
        eb_coeff[e_below_key] += 1.0
    else # if this edge points away from the center (+r direction, same as the current edge)
        eb_coeff[e_below_key] += -1.0
    end

    return eb_coeff
end


"""
    a_angle_coeff(primalmesh::Primalmeshstruct, eb_current::Array{Any,1}, fbound2::Array{Any,1}, ebound2::Array{Any,1},
    r2::Float64, r3::Float64, Lmax::Int64)

    calculate the coefficients to be put in the matrix for open BC of the field tangential to the boundary sphere
"""
function a_angle_coeff(primalmesh::Primalmeshstruct, eb_current::SArray{Tuple{2},Int64,1,2}, fbound2::Array{Any,1}, ebound2::Array{Any,1},
                       r2::Float64, r3::Float64, Lmax::Int64, Lmin::Int64, k::Complex{Float64})

    # info about the vertices of this edge
    n1  = primalmesh.nodedict[eb_current[1]].coords
    n2  = primalmesh.nodedict[eb_current[2]].coords
    _, theta_v1, phi_v1 = findSphericalCoords(n1[1],n1[2],n1[3])
    _, theta_v2, phi_v2 = findSphericalCoords(n2[1],n2[2],n2[3])
    r_unit_v1 = [sin(theta_v1)*cos(phi_v1), sin(theta_v1)*sin(phi_v1), cos(theta_v1)]
    r_unit_v2 = [sin(theta_v2)*cos(phi_v2), sin(theta_v2)*sin(phi_v2), cos(theta_v2)]

    vec12 =  n2-n1
    # find the vector spherical Harmonics at the two vertices
    Psi_lm_eb = Dict{SVector{2, Int}, Complex{Float64}}()
    Phi_lm_eb = Dict{SVector{2, Int}, Complex{Float64}}()
    for l = Lmin:Lmax
        for m = -l:l
            Y_lm_n1 = computeYlm(theta_v1, phi_v1, lmax = l)[(l,m)]
            Y_lm_n2 = computeYlm(theta_v2, phi_v2, lmax = l)[(l,m)]
            if m<l
                Y_lmp1_n1 = computeYlm(theta_v1, phi_v1, lmax = l)[(l,m+1)] 
                Y_lmp1_n2 = computeYlm(theta_v2, phi_v2, lmax = l)[(l,m+1)] 
            else
                Y_lmp1_n1 = 0.0
                Y_lmp1_n2 = 0.0
            end
            # this is actually Psi_lm/r (= grad_Y_lm) since r is absorbed into the expansion coefficients of the field
            Psi_lm_1 = (1/r2)*[cos(theta_v1)*cos(phi_v1)*(m*cot(theta_v1)*Y_lm_n1 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v1)*Y_lmp1_n1) - im*sin(phi_v1)/sin(theta_v1)*m*Y_lm_n1,
                               cos(theta_v1)*sin(phi_v1)*(m*cot(theta_v1)*Y_lm_n1 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v1)*Y_lmp1_n1) + im*cos(phi_v1)/sin(theta_v1)*m*Y_lm_n1,
                               -sin(theta_v1)*(m*cot(theta_v1)*Y_lm_n1 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v1)*Y_lmp1_n1)]
            Psi_lm_2 = (1/r2)*[cos(theta_v2)*cos(phi_v2)*(m*cot(theta_v2)*Y_lm_n2 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v2)*Y_lmp1_n2) - im*sin(phi_v2)/sin(theta_v2)*m*Y_lm_n2,
                               cos(theta_v2)*sin(phi_v2)*(m*cot(theta_v2)*Y_lm_n2 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v2)*Y_lmp1_n2) + im*cos(phi_v2)/sin(theta_v2)*m*Y_lm_n2,
                               -sin(theta_v2)*(m*cot(theta_v2)*Y_lm_n2 + sqrt((l-m)*(l+m+1))*exp(-im*phi_v2)*Y_lmp1_n2)]
            # this is actually Phi_lm/r (= unit_r x grad_Y_lm) since r is absorbed into the expansion coefficients of the field
            Phi_lm_1 = cross(r_unit_v1, Psi_lm_1)
            Phi_lm_2 = cross(r_unit_v2, Psi_lm_2)
            # Psi_lm and Phi_lm values along the edge (1-form, not vector)
            Psi_lm_eb[[l,m]] = 0.5*dot(Psi_lm_1+Psi_lm_2, vec12)
            Phi_lm_eb[[l,m]] = 0.5*dot(Phi_lm_1+Phi_lm_2, vec12)
        end
    end

    # prepping the list of hankel SpecialFunctions
    hankelH1_1_2  = Dict{Int, Complex{Float64}}() # J_1/2
    hankelH1_m1_2 = Dict{Int, Complex{Float64}}() # J_-1/2
    hankelH1_3_2  = Dict{Int, Complex{Float64}}() # J_3/2
    hankelH1_m3_2 = Dict{Int, Complex{Float64}}() # J_-3/2
    hankelH1_5_2  = Dict{Int, Complex{Float64}}() # J_5/2
    for l = Lmin:Lmax
        hankelH1_1_2[l]  = hankelh1(l+1/2,k*r2)
        hankelH1_m1_2[l] = hankelh1(l-1/2,k*r2)
        hankelH1_3_2[l]  = hankelh1(l+3/2,k*r2)
        hankelH1_m3_2[l] = hankelh1(l-3/2,k*r2)
        hankelH1_5_2[l]  = hankelh1(l+5/2,k*r2)
    end

    # the dictionary to store the coefficients (the keys are the edges that contribute to the boundary condition)
    ###eb_coeff = Dict{SVector{2, Int}, Float64}()
    eb_coeff = Dict{Array{Any,1}, Complex{Float64}}()
    dr23 = r3-r2

    for ekey in ebound2
        # find the two boundary faces that share this edge
        facelist  = fbound2[findall(x->(ekey[1] in x)&&(ekey[2] in x), fbound2)]
        # info on vertices of the ekey
        v1_coords = primalmesh.nodedict[ekey[1]].coords
        v2_coords = primalmesh.nodedict[ekey[2]].coords
        edge_linepoint = [v1_coords, v2_coords]
        f_contri_to_eb = 0.0
        for fkey in facelist
            temp_dvec = zeros(3)
            area = primalmesh.facedict[fkey].area
            circ = get_circumcenter_face(primalmesh.facedict[fkey], primalmesh.nodedict)
            edgelist   = primalmesh.facedict[fkey].edges
            tri_verts  = [primalmesh.nodedict[fkey[1]].coords, primalmesh.nodedict[fkey[2]].coords, primalmesh.nodedict[fkey[3]].coords]
            ###tri_angles = triAngleCompute(tri_verts[1], tri_verts[2], tri_verts[3])
            # check whether the circumcenter is inside the triangle 
            relpos = insideTriCheck(circ, tri_verts) # = true or false
            # find the remaining vertex of the triangle fkey that is not a vertex of ekey
            v_rem = fkey[findfirst(x->!(x in ekey), fkey)]
            # check whether the circumcenter and the remaining vertex lie in the same half space created by edge ekey
            reldir_v_rem = halfspace_line_check(circ, primalmesh.nodedict[v_rem].coords, edge_linepoint) # = true or false
            # find the vector normal to the current edge ekey and points towards the opposite vertex (v_rem)
            edge_center   = 0.5*(v1_coords + v2_coords)
            edge2circ     = (-1+2*reldir_v_rem)*(circ - edge_center) 
            edge2circ_len = norm(edge2circ)
            # weighted contribution from each node of the edge ekey
            for vert in ekey
                # the other vertex of the edge
                v_other = ekey[findfirst(x->x!=vert,ekey)]
                # find the edge opposite to v_other
                op_edge = edgelist[findfirst(x->!(v_other in x), edgelist)]
                # height from vertices to opposite edges
                h = 2.0*area/primalmesh.edgedict[op_edge].length
                # check whether the circumcenter and v_other lie in the same half space created by edge op_edge
                op_edge_linepoint = [primalmesh.nodedict[op_edge[1]].coords, primalmesh.nodedict[op_edge[2]].coords]
                reldir_v_other = halfspace_line_check(circ, primalmesh.nodedict[v_other].coords, op_edge_linepoint) # = true or false
                # find the vector normal to the opposite edge and points towards the vertex v_other
                op_edge_center   = 0.5*(primalmesh.nodedict[op_edge[1]].coords + primalmesh.nodedict[op_edge[2]].coords) 
                op_edge2circ     = (-1+2*reldir_v_other)*(circ - op_edge_center) 
                op_edge2circ_len = norm(op_edge2circ)     
                op_edge_norm     = op_edge2circ/op_edge2circ_len # normalized
                # compute the weight based on the subdivision of the triangle
                if relpos 
                    w  = 0.25*(primalmesh.edgedict[op_edge].length*op_edge2circ_len + primalmesh.edgedict[ekey].length*edge2circ_len)/area
                else # if circ not inside triangle, have to calculate intersections explicitly
                    # edge that doesn't include vert (the edge opposite to vert)
                    novert_edge      = edgelist[findfirst(x->!(vert in x), edgelist)]
                    linepoint_novert = [primalmesh.nodedict[novert_edge[1]].coords, primalmesh.nodedict[novert_edge[2]].coords] 
                    # check whether vert and the circumcenter are in the same halfspace created by this edge 
                    reldir_vert = halfspace_line_check(circ, primalmesh.nodedict[vert].coords, linepoint_novert) # = true or false
                    if !(reldir_vert)
                        # distance form circ to the edge opposite to vert
                        novert_edge2circ_len = norm(circ - 0.5*(primalmesh.nodedict[v_other].coords+primalmesh.nodedict[v_rem].coords))
                        v_other_angle, v_rem_angle, vert_angle = triAngleCompute(Array{Float64,1}(primalmesh.nodedict[v_other].coords), Array{Float64,1}(primalmesh.nodedict[v_rem].coords), 
                                                                                 Array{Float64,1}(primalmesh.nodedict[vert].coords))
                        seg1 = 0.5*primalmesh.edgedict[ekey].length/cos(v_other_angle)
                        seg2 = 0.5*primalmesh.edgedict[op_edge].length/cos(v_rem_angle)
                        appendix_area = 0.5*(novert_edge2circ_len*(primalmesh.edgedict[novert_edge].length-seg1-seg2))
                        w = (0.25*(primalmesh.edgedict[op_edge].length*op_edge2circ_len + primalmesh.edgedict[ekey].length*edge2circ_len) - appendix_area)/area
                    else
		        v_other_angle, v_rem_angle, vert_angle = triAngleCompute(Array{Float64,1}(primalmesh.nodedict[v_other].coords), Array{Float64,1}(primalmesh.nodedict[v_rem].coords), 
                                                                                 Array{Float64,1}(primalmesh.nodedict[vert].coords))
                        # check whether v_rem and circ are in the same halfspace created by ekey
                        if reldir_v_rem 
                            w = 0.5*(0.5*primalmesh.edgedict[ekey].length)^2*tan(vert_angle)/area
                        else
                            w = 0.5*(0.5*primalmesh.edgedict[op_edge].length)^2*tan(vert_angle)/area
                        end
                    end
                end
                # find the orientation of the edge ekey
                if vert == ekey[1]
                    sgn = 1.0
                else
                    sgn = -1.0
                end
                temp_dvec += w*sgn*(1/h)*op_edge_norm
            end
            temp_dvec = (area/r2^2)*temp_dvec

            _, theta_f, phi_f = findSphericalCoords(circ[1], circ[2], circ[3])
            r_unit_f = [sin(theta_f)*cos(phi_f), sin(theta_f)*sin(phi_f), cos(theta_f)]
            r_lm_chunk = [0.0, 0.0, 0.0]
            for l = Lmin:Lmax
                for m = -l:l
                    Y_lm = computeYlm(theta_f, phi_f, lmax = l)[(l,m)]
                    if m<l
                        Y_lmp1 = computeYlm(theta_f, phi_f, lmax = l)[(l,m+1)]  
                    else
                        Y_lmp1 = 0.0
                    end
                    grad_Ylm_f = (1/r2)*[cos(theta_f)*cos(phi_f)*(m*cot(theta_f)*Y_lm + sqrt((l-m)*(l+m+1))*exp(-im*phi_f)*Y_lmp1) - im*sin(phi_f)/sin(theta_f)*m*Y_lm,
                                        cos(theta_f)*sin(phi_f)*(m*cot(theta_f)*Y_lm + sqrt((l-m)*(l+m+1))*exp(-im*phi_f)*Y_lmp1) + im*cos(phi_f)/sin(theta_f)*m*Y_lm,
                                        -sin(theta_f)*(m*cot(theta_f)*Y_lm + sqrt((l-m)*(l+m+1))*exp(-im*phi_f)*Y_lmp1)]
                    Psi_lm_f = r2*grad_Ylm_f
                    Phi_lm_f = r2*cross(r_unit_f, grad_Ylm_f)

                    r_lm_chunk += 1/(l*(l+1))*(
                                               conj(Phi_lm_f)*(r2 + dr23*(0.5 + k*r2*0.5*(hankelH1_m1_2[l]-hankelH1_3_2[l])/hankelH1_1_2[l]))*Phi_lm_eb[[l,m]]
                                              +conj(Psi_lm_f)*(r2 + dr23*(1.0 + 0.5*(-3.0*hankelH1_1_2[l] + (k*r2)^2*(hankelH1_m3_2[l] - 2.0*hankelH1_1_2[l] + hankelH1_5_2[l]))/(hankelH1_1_2[l] + k*r2*(hankelH1_m1_2[l]-hankelH1_3_2[l]))))*Psi_lm_eb[[l,m]] 
                                              )
                end
            end
            f_contri_to_eb += dot(temp_dvec, r_lm_chunk)
        end
        eb_coeff[ekey] = f_contri_to_eb
    end

    return eb_coeff
end
