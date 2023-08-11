"""
    doublecurl(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound::Array{Any,1})

    compute curl-curl operator
"""
function doublecurl(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound::Array{Any,1})
    Ne = length(primalmesh.edgedict)
    dcurl = Array{Float64,2}(undef, Ne, Ne)
    
    edgekeylist = collect(keys(primalmesh.edgedict))
    facekeylist = collect(keys(primalmesh.facedict))

    for i = 1:Ne
        ki = edgekeylist[i]
        if !(ki in ebound)
            v1 = ki[1]
            v2 = ki[2]
            dual_area  = dualmesh.dualfacedicts.interior_dualfacedict[ki].raw_area
            primal_len = primalmesh.edgedict[ki].length
            for j = 1:Ne
                kj = edgekeylist[j]
                if (v1==kj[1])&&(v2!=kj[2]) # if the current edge and neighboring edge share v1 node
                    facekey     = sort([v1,v2,kj[2]])
                    if (facekey in facekeylist)
                        area        = primalmesh.facedict[facekey].area
                        dual_len    = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j]  = -(dual_len/area)*(primal_len/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop       
                        dcurl[i,i] += (dual_len/area)*(primal_len/dual_area)  # for each loop coming out of v1, add up contribution of the central edge
                    end
                elseif (v1==kj[2])&&(v2!=kj[1]) # if the current edge and neighboring edge share v2 node
                    facekey    = sort([v1,v2,kj[1]])
                    if (facekey in facekeylist)
                        area        = primalmesh.facedict[facekey].area
                        dual_len    = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j]  = (dual_len/area)*(primal_len/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                        dcurl[i,i] += (dual_len/area)*(primal_len/dual_area) # for each loop coming out of v1, add up contribution of the central edge
                    end
                elseif (v2==kj[1])&&(v1!=kj[2])
                    facekey    = sort([v1,v2,kj[2]])
                    if (facekey in facekeylist)
                        area       = primalmesh.facedict[facekey].area
                        dual_len   = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j] = (dual_len/area)*(primal_len/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                    end
                elseif (v2==kj[2])&&(v1!=kj[1])
                    facekey    = sort([v1,v2,kj[1]])
                    if (facekey in facekeylist)
                        area       = primalmesh.facedict[facekey].area
                        dual_len   = dualmesh.dualedgedicts.interior_dualedgedict[facekey].length
                        dcurl[i,j] = -(dual_len/area)*(primal_len/dual_area) # the sign (+ or -) needs to be taken cared of to close the loop
                    end
                end
            end
        end
    end
    
    return dcurl
end


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
