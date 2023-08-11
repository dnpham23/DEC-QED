include("../../../gmesh_parser/dec-qed/julia_program/dualmesh.jl");
include("utils.jl");
include("mappings.jl");
using JLD2;

primalmesh = complete_primalmesh(raw"/meshes/single_cube_dense.msh");
dualmesh   = complete_dualmesh(raw"/meshes/single_cube_dense.msh");
Ne         = length(primalmesh.edgedict);
Nv         = length(primalmesh.nodedict);

@load "SimulationData3.jld2" phi Eigs xmin xmax ymin ymax zmin zmax lambda;

Neig = length(Eigs);
phi_vec = zeros(Nv, 3, Neig); 

for n = 1:Neig
    phi_vec[:,:,n] = pp_sharp(primalmesh, dualmesh, real(phi[:,n]), Nv, Ne) 

end



@save "SimulationData.jld2" Eigs phi xmin xmax ymin ymax zmin zmax lambda phi_vec
