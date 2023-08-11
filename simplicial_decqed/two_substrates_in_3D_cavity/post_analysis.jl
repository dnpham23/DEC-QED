include("../../../../gmesh_parser/dec-qed/julia_program/dualmesh.jl");
include("utils.jl");
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

# find all edges that live on the 2d surface (that mounts the qubit) of substrate 1
substr1_surface_edgelist = find_substr1_surface_edge(primalmesh.edgedict, physicalnames_dict, all_entities_struct);
# still missing in substr1_surface_edgelist: edges that has just one vertex on the JJ
# first, find the nodes on JJ1
JJ1_nodelist = find_JJ1_nodes(primalmesh.nodedict, physicalnames_dict, all_entities_struct);
# collect the remaining surface edges attached to these nodes 
JJ_attached_edges = [];
for e in collect(keys(primalmesh.edgedict))
    if e[1] in JJ1_nodelist 
        if (abs(primalmesh.nodedict[e[2]].coords[3] - t/2)/(t/2) <= tol) # if this edge lies on substrate 1's surface
            if !(e in substr1_surface_edgelist)
                push!(substr1_surface_edgelist, e)
                push!(JJ_attached_edges, e)
            end
        end
    elseif e[2] in JJ1_nodelist
        if (abs(primalmesh.nodedict[e[1]].coords[3] - t/2)/(t/2) <= tol) # if this edge lies on substrate 1's surface
            if !(e in substr1_surface_edgelist)
                push!(substr1_surface_edgelist, e)
                push!(JJ_attached_edges, e)
            end
        end
    end
end

# find all the coordinates of the endpoints of these edges
substr1_surface_edge_endpt_coords = [];
for e in substr1_surface_edgelist
    n1 = e[1]
    n2 = e[2]
    push!(substr1_surface_edge_endpt_coords, [collect(primalmesh.nodedict[n1].coords), collect(primalmesh.nodedict[n2].coords)])
end

# find all the coordinates of the endpoints of surface edges pointing out of the JJ
JJ_attached_edges_endpt_coords = [];
for e in JJ_attached_edges
    n1 = e[1]
    n2 = e[2]
    push!(JJ_attached_edges_endpt_coords, [collect(primalmesh.nodedict[n1].coords), collect(primalmesh.nodedict[n2].coords)])
end


# load data for plotting
@load("SimulationData3.jld2");
mode_id = 10; # index of mode to be plotted
# collect the fields on substrate 1's surface only
substr1_surface_edge_ind = findall(x->x in substr1_surface_edgelist, edgekeylist);
phi_substrate1_surface = phi[substr1_surface_edge_ind,:];
# normalize the edge field
maxphi_mode_id = findmax(abs.(phi_substrate1_surface[:,mode_id]))[1];
phi_substrate1_surface_mode_id_normalized = real.(phi_substrate1_surface[:,mode_id])./maxphi_mode_id;
phi_substrate1_surface_mode_id_normalized = abs.(phi_substrate1_surface_mode_id_normalized); # eliminate the signs due to orientations of the edges
# 
base = 0.15;
phi_hot = phi_substrate1_surface_mode_id_normalized[findall(x->x>base, phi_substrate1_surface_mode_id_normalized)];
phi_hot = phi_hot .+ (1.0 .- phi_hot).^3.5;
phi_substrate1_surface_mode_id_normalized[findall(x->x>base, phi_substrate1_surface_mode_id_normalized)] = phi_hot;

#phi_substrate1_surface_mode_id_normalized = (phi_substrate1_surface_mode_id_normalized.- 0.5).*2;

using PyPlot
# plot the edge fields
fig = figure(figsize=(c1/2, (d+2*c2)/2))
for i in 1:length(substr1_surface_edgelist)
    coord1 = substr1_surface_edge_endpt_coords[i][1]
    coord2 = substr1_surface_edge_endpt_coords[i][2]
    f = phi_substrate1_surface_mode_id_normalized[i]
    if f > base
        plot([coord1[1], coord2[1]], [coord1[2], coord2[2]], color=(f, 0, 0),linewidth = 2.0)
    else
        plot([coord1[1], coord2[1]], [coord1[2], coord2[2]], color=(0, 0, 1-f*4),linewidth = 2.0)
    end
end
xlabel("x"); ylabel("y");


# plot the fields on edges close to JJ
JJ_attached_edges_ind = findall(x->x in JJ_attached_edges, edgekeylist);
for i = 1:length(JJ_attached_edges)
    coord1 = JJ_attached_edges_endpt_coords[i][1]
    coord2 = JJ_attached_edges_endpt_coords[i][2]
    plot([coord1[1], coord2[1]], [coord1[2], coord2[2]], color=(0.75 + (rand(1)[1].- 0.5)/2, 0, 0),linewidth = 2.0)
end