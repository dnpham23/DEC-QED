using Distributed

# launch worker processes
num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(num_cores)


@everywhere begin
    include("../../../../gmesh_parser/dec-qed/julia_program/dualmesh.jl");
    include("sphericalmesh.jl");
    include("utils.jl");
    include("mappings.jl");
    include("openbc.jl");
    include("operators.jl");
    using SharedArrays;
end

# obtain information about primal mesh and construct the dual mesh
#primalmesh, physicalnames_dict, all_entities_struct = complete_primalmesh(raw"/hankel/meshes/cocentric_spheres_dielectric_sphere_2.msh");
#dualmesh, primalmesh   = complete_dualmesh(raw"/hankel/meshes/cocentric_spheres_dielectric_sphere_lessdense.msh");
primalmesh, physicalnames_dict, all_entities_struct = complete_primalmesh(raw"/hankel/meshes/cocentric_spheres_dielectric_sphere_lessdense.msh");
dualmesh, primalmesh   = complete_dualmesh(raw"/hankel/meshes/cocentric_spheres_dielectric_sphere_lessdense.msh");
# numerical error tolerance for coordinates of vertices
tol = 1e-3;
# the three spherical boundary layers and their origin 
r1  = 16.0;
r2  = 20.0;
r3  = 24.0;
rd  = 12.0;
xo  = 0.0;
yo  = 0.0;
zo  = 0.0;
# basis size 
Lmax = 10; # maximum order of bessel functions
Lmin = 1;
# dielectric constant
dconst = 3.0;

# find boundary nodes and edges, then add another layer
vbound, vbound1, vbound2, vbound3, Nvbound_all = vboundary_sphere(r1, r2, r3, xo, yo, zo, primalmesh, tol);
ebound1, ebound2, ebound3, ebound_r12, ebound_r23 = eboundary_sphere(vbound1, vbound2, vbound3, primalmesh, r1, r2, r3, xo, yo, zo, tol);
fbound1, fbound2, fbound3 = fboundary_sphere(vbound1, vbound2, vbound3, primalmesh);
nodedict_extra, vbound4, r4 = bnode_add(vbound3, primalmesh, r1, r2, r3, xo, yo, zo, tol);
edgedict_extra = bedge_add(vbound3, vbound4, primalmesh, nodedict_extra, r3, r4);
ebound_extra   = collect(keys(edgedict_extra));

Ne          = length(primalmesh.edgedict) + length(edgedict_extra);
Nv          = length(primalmesh.nodedict) + length(nodedict_extra);
edgekeylist = [collect(keys(primalmesh.edgedict)); ebound_extra]; # newly contructed boundary edges included
nodekeylist = [collect(keys(primalmesh.nodedict)); vbound4]; # newly contructed boundary nodes included
ebound_ind  = findall(x->x in [ebound3; ebound_extra], edgekeylist);
edge_dielectric, edge_dielec_ind = dielectric_edgecollect(primalmesh, rd, tol, xo, yo, zo, edgekeylist);

# the operators 
dcurl = doublecurl_openbc(primalmesh, dualmesh, ebound3, ebound_extra, edgedict_extra, edgekeylist);

# range of wave numbers to sweep over
k_re = collect(LinRange(-11.875*pi/r2, -0.125*pi/r2, 48));
k_im = collect(LinRange(-11.875*pi/r2, -0.125*pi/r2, 48));
klist = [];
for kr in k_re
    for ki in k_im
        append!(klist, kr+im*ki)
    end
end

Nsv = 15;
svlist = SharedArray{Float64}((length(klist), Nsv));

# boundary operator
@sync @distributed for i = 1:length(klist)
    k = klist[i]
    println("step = ", i)
    diag_array = vec(ones(Ne))
    diag_array[edge_dielec_ind] .= dconst^2; # modify the wavenumber in the dielectric region
    diag_array[ebound_ind]   .= 0.0; # filter out edges at the boundary
    Ik2 = k^2*diagm(diag_array) 
    bop = boundary_op(primalmesh, fbound2, ebound2, ebound3, ebound_r12, ebound_r23, ebound_extra, vbound2, r1, r2, r3, Lmax, Lmin, k, edgekeylist, Ne)
    mop = dcurl + bop - Ik2
    svlist[i,:] = svdvals(mop)[end-Nsv+1:end]
end

meshfile = "cocentric_spheres_dielectric_sphere_lessdense.msh";

using JLD2
@save "SimulationData.jld2" klist svlist meshfile Lmax Lmin Nsv r1 r2 r3 rd dconst
###################################
