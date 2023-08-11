include("../../../../gmesh_parser/dec-qed/julia_program/dualmesh.jl");
include("utils.jl");
include("operators.jl");
include("mappings.jl");

# obtain information about primal mesh and construct the dual mesh
primalmesh, physicalnames_dict, all_entities_struct = complete_primalmesh(raw"/meshes/two_JJs_on_substr_in_3d_cav.msh");
dualmesh, primalmesh   = complete_dualmesh(raw"/meshes/two_JJs_on_substr_in_3d_cav.msh");
Ne         = length(primalmesh.edgedict);
Nv         = length(primalmesh.nodedict);

### geometrical properties
# Parameters for the 3d cavity
L = 1000;
H = 550;
R = 80;
# Parameters for the two identical substrates
s = 50;
w = 25;
t = 5;
l = 600;
# Parameters for the JJ pads
c1 = 12;
c2 = 10;
d  = 8;
# material properties
λₗ = 0.1; # London penetration depth of the superconductor
ϵ = 11.68; # dielectric constant of the substrate

# numerical error tolerance for coordinates of vertices
tol = 1e-5;

# find all edges that belongs to the superconducting devices
sc_edgelist = find_sc_edges(primalmesh.edgedict, primalmesh.nodedict, physicalnames_dict, all_entities_struct);

# find all edges that belongs to the substrate
substr_edgelist = find_substr_edges(primalmesh.edgedict, physicalnames_dict, all_entities_struct, sc_edgelist);



# find dimensions of the cavity 
allcoords = getfield.(getindex.(Ref(primalmesh.nodedict), keys(primalmesh.nodedict)), :coords);
xmin      = findmin([x[1] for x in allcoords])[1];
xmax      = findmax([x[1] for x in allcoords])[1];
ymin      = findmin([x[2] for x in allcoords])[1];
ymax      = findmax([x[2] for x in allcoords])[1];
zmin      = findmin([x[3] for x in allcoords])[1];
zmax      = findmax([x[3] for x in allcoords])[1];

# find boundary edges
ebound = eboundary_gm(primalmesh.edgedict, physicalnames_dict, all_entities_struct);

# compute differential operators
dcurl = doublecurl(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound::Array{Any,1});

# list of edge keys
edgekeylist = collect(keys(primalmesh.edgedict));
nodekeylist = collect(keys(primalmesh.nodedict));

# apply material properties to the edges 
invlambda2 = zeros(Ne);
invepsilon =  ones(Ne);
invlambda2[findall(x->x in sc_edgelist, edgekeylist)]  .= 1/λₗ^2;
invepsilon[findall(x->x in substr_edgelist, edgekeylist)] .= 1/ϵ;

ebound_ind  = findall(x->x in ebound, edgekeylist);
invlambda2[ebound_ind] .= 0.0; # filter out edges at the boundary
invepsilon[ebound_ind] .= 0.0; # filter out edges at the boundary

# final matrix
A = dcurl.*invepsilon + Diagonal(invlambda2);
###A = dcurl;

# apply Dirichlet boundary condition
for ind in ebound_ind
    A[ind,ind] = 1.0e5;
end

using LinearAlgebra
using Arpack


Neig      = 10;
sigma     = 1.4*pi^2;
IMat      = Diagonal(vec(ones(Ne)));
Eigs, phi = eigs(A, IMat, nev = Neig, which=:LM, maxiter=1000, sigma = sigma);
###Eigs, phi = eigs(A, nev = Neig, which=:SM, sigma = pi^2);
###Eigs, phi = eigs(A, nev = Neig, which=:LM, maxiter=1000, sigma = 3*pi^2);


# compute the vector field
phi_vec = zeros(Nv, 3, Neig); 

for n = 1:Neig
    phi_vec[:,:,n] = pp_sharp(primalmesh, dualmesh, real(phi[:,n]), Nv, Ne) 

end

nodecoords = zeros(Nv, 3);
for nv = 1:Nv
    nodecoords[nv,:] = primalmesh.nodedict[nodekeylist[nv]].coords
end

meshfile = "two_JJs_on_substr_in_3d_cav.msh"

using JLD2;
@save "SimulationData3.jld2" Eigs phi xmin xmax ymin ymax zmin zmax t d λₗ ϵ phi_vec nodecoords meshfile sigma 


