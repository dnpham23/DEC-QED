using JLD2
@load("W100_Tm1000_u0_55_alpha0_beta0.jld2");
using Plots
gr(size = (1500, 600), legend = false, colorbar = false)

endcut = nv_t;


phimax = findmax(phi)[1];
phimin = findmin(phi)[1];
anim = @animate for i = 1:36:endcut
                  plot(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phi[i,:], linewidth = 2, xlabel = "x", ylabel="\$\\phi\$", ylims = (phimin*1.1, phimax*1.1));
       end
gif(anim, "phi_vs_time.gif", fps = 40)

phit = phi[2:end,:] .- phi[1:end-1,:];
phitmax = findmax(phit)[1];
phitmin = findmin(phit)[1];
anim = @animate for i = 1:36:endcut-1
                  plot(-(nv_x-1)/2*dx:dx:(nv_x-1)/2*dx, phit[i,:], linewidth = 2, xlabel = "x", ylabel="\$\\phi_t\$", ylims = (1.1*phitmin, 1.1*phitmax));
              end
phit = 0;
gif(anim, "phit_vs_t.gif", fps = 40)

phix = phi[:,2:end] .- phi[:,1:end-1];
phixmax = findmax(phix)[1];
phixmin = findmin(phix)[1];
anim = @animate for i = 1:36:endcut
                  plot(-(nv_x-1)/2*dx+dx/2:dx:(nv_x-1)/2*dx-dx/2, phix[i,:], linewidth = 2, xlabel = "x", ylabel="\$\\phi_x\$", ylims = (1.1*phixmin, 1.1*phixmax));
              end
phix = 0;
gif(anim, "phix_vs_t.gif", fps = 40)





