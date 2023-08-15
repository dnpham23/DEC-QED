include("../gmsh_parser/dualmesh.jl");
include("utils.jl");
include("mappings.jl");
include("operators.jl");

# obtain information about primal mesh and construct the dual mesh
primalmesh, physicalnames_dict, all_entities_struct = complete_primalmesh(raw"/meshes/embedded_rectangle_dense.msh");
dualmesh, primalmesh   = complete_dualmesh(raw"/meshes/embedded_rectangle_dense.msh");
Ne         = length(primalmesh.edgedict);
Nv         = length(primalmesh.nodedict);

# list containing the penetration depths of the materials in the system, in the same order as how volumes are defined in .msh file
λₗ = [1e4, 0.3] 
invlambda2 = 1 ./λₗ.^2

# find effective value of 1/λₗ² at the primal edges (finding effective areas of dual faces)
invlambda2_eff = effective_material(primalmesh.nodedict, primalmesh.edgedict, primalmesh.facedict, primalmesh.tetdict, invlambda2, dualmesh.dualfacedicts);

# numerical error tolerance for coordinates of vertices
tol = 1e-5;

# find dimensions of the cavity 
allcoords = getfield.(getindex.(Ref(primalmesh.nodedict), keys(primalmesh.nodedict)), :coords);
xmin      = findmin([x[1] for x in allcoords])[1];
xmax      = findmax([x[1] for x in allcoords])[1];
ymin      = findmin([x[2] for x in allcoords])[1];
ymax      = findmax([x[2] for x in allcoords])[1];
zmin      = findmin([x[3] for x in allcoords])[1];
zmax      = findmax([x[3] for x in allcoords])[1];

# find boundary nodes
vbound = vboundary(xmax,xmin,ymax,ymin,zmax,zmin, primalmesh, tol);
ebound = eboundary(vbound, primalmesh, tol);

# compute differential operators
dcurl = doublecurl(primalmesh::Primalmeshstruct, dualmesh::Dualmeshstruct, ebound::Array{Any,1});

# material property
edgekeylist = collect(keys(primalmesh.edgedict));
ebound_ind  = findall(x->x in ebound, edgekeylist);
invlambda2_eff[ebound_ind] .= 0.0; # filter out edges at the boundary

# final matrix
A = dcurl + Diagonal(invlambda2_eff);
###A = dcurl;

# apply Dirichlet boundary condition
for ind in ebound_ind
    A[ind,ind] = 1.0e5;
end

using LinearAlgebra
using Arpack

Neig      = 10;
IMat      = Diagonal(vec(ones(Ne)));
Eigs, phi = eigs(A, IMat, nev = Neig, which=:LM, maxiter=1000, sigma = 5*pi^2);
###Eigs, phi = eigs(A, nev = Neig, which=:SM, sigma = pi^2);
###Eigs, phi = eigs(A, nev = Neig, which=:LM, maxiter=1000, sigma = 3*pi^2);

# compute the vector field
phi_vec = zeros(Nv, 3, Neig); 

for n = 1:Neig
    phi_vec[:,:,n] = pp_sharp(primalmesh, dualmesh, real(phi[:,n]), Nv, Ne) 

end



using JLD2;
@save "SimulationData3.jld2" Eigs phi xmin xmax ymin ymax zmin zmax lambda phi_vec 