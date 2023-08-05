function ampere_timedomain(phixpast1::Array{Float64,1},phixpast2::Array{Float64,1},phiypast1::Array{Float64,1},phiypast2::Array{Float64,1}, psixnow::Array{Float64,2},psiynow::Array{Float64,2},psixpast1::Array{Float64,1},psiypast1::Array{Float64,1}, ecurrent::Array{Float64,2}, ebound::Array{Int64,1}, dt::Float64, Ne_airx::Int64, Ne_airy::Int64, Ne_x::Int64, lx::Float64, ly::Float64, eps::Float64, mu::Float64,A1::Float64,A2::Float64,e_air::Array{Int,1},e_sc::Array{Int,1},ne_airx::Int64,ne_airy::Int64)
    
    phix = zeros(Ne_airx,1);
    phiy = zeros(Ne_airy,1);
    Nebound_sc  = length(ebound) + ne_airx*2 + ne_airy*2;
    ebound_sc   = [ebound; e_sc[end-(ne_airx*2+ne_airy*2)+1:end]];
    e_interface = e_sc[end-(ne_airx*2+ne_airy*2)+1:end];
    for i = 1:Ne_airx # index in the air array only (not global index)
        real_ind = e_air[i]; # actual global edge index
        if !(real_ind in ebound_sc)
            rownum  = div(i-1,ne_airx) + 1;
            colnum  = (i-1)%ne_airx + 1;
            yedge   = (rownum-2)*(ne_airx+1) + colnum; # index in the array for y-edge in air (not global edge index)
            phix[i] = A2*(lx/ly)*(dt^2)/eps*ecurrent[real_ind] - (psixnow[real_ind]-psixpast1[real_ind])*dt - A1*(dt^2)/(mu*eps*ly^2)*(-phixpast1[i-ne_airx] + 2*phixpast1[i] -phixpast1[i+ne_airx] + phiypast1[yedge] - phiypast1[yedge+1] - phiypast1[yedge+ne_airx+1] + phiypast1[yedge+ne_airx+2]) + 2*phixpast1[i] - phixpast2[i];
            
        elseif (real_ind in e_interface) # if the edge is at the interface
            rownum  = div(i-1,ne_airx) + 1;
            colnum  = (i-1)%ne_airx + 1;
            yedge   = (rownum-2)*(ne_airx+1) + colnum; # index in the array for y-edge in air (not global edge index)
            if rownum == 1 # if the edge is at the lower horizontal interface
                phix[i] = A2*(lx/ly)*(dt^2)/eps*ecurrent[real_ind] - (psixnow[real_ind]-psixpast1[real_ind])*dt - A1*(dt^2)/(mu*eps*ly^2)*( phixpast1[i] -phixpast1[i+ne_airx] - phiypast1[yedge+ne_airx+1] + phiypast1[yedge+ne_airx+2]) + 2*phixpast1[i] - phixpast2[i];
            elseif (rownum == ne_airy+1) # if the edge is at the upper horizontal interface
                phix[i] = A2*(lx/ly)*(dt^2)/eps*ecurrent[real_ind] - (psixnow[real_ind]-psixpast1[real_ind])*dt - A1*(dt^2)/(mu*eps*ly^2)*(-phixpast1[i-ne_airx] + phixpast1[i] + phiypast1[yedge] - phiypast1[yedge+1]) + 2*phixpast1[i] - phixpast2[i];
            end
        end
    end
    
    for i = 1:Ne_airy # index in the air array only (not global index)
        real_ind = e_air[i+Ne_airx]; # actual global edge index
        if !(real_ind in ebound_sc)
            rownum = div(i-1,ne_airx+1) + 1;
            colnum = (i-1)%(ne_airx+1) + 1;
            xedge = ne_airx*(rownum-1) + colnum-1; # index in the array for x-edge in air (not global edge index)
            phiy[i] = A2*(ly/lx)*(dt^2)/eps*ecurrent[real_ind] - (psiynow[real_ind-Ne_x]-psiypast1[real_ind-Ne_x])*dt - A1*(dt^2)/(mu*eps*lx^2)*(-phiypast1[i-1] + 2*phiypast1[i] - phiypast1[i+1] + phixpast1[xedge] - phixpast1[xedge+ne_airx] - phixpast1[xedge+1] + phixpast1[xedge+ne_airx+1]) + 2*phiypast1[i] - phiypast2[i];
                
        elseif (real_ind in e_interface) # if the edge is at the interface
            rownum = div(i-1,ne_airx+1) + 1;
            colnum = (i-1)%(ne_airx+1) + 1;
            xedge = ne_airx*(rownum-1) + colnum-1; # index in the array for x-edge in air (not global edge index)
            if colnum ==1  # if the edge is at the left vertical interface
                phiy[i] = A2*(ly/lx)*(dt^2)/eps*ecurrent[real_ind] - (psiynow[real_ind-Ne_x]-psiypast1[real_ind-Ne_x])*dt - A1*(dt^2)/(mu*eps*lx^2)*(phiypast1[i] - phiypast1[i+1] - phixpast1[xedge+1] + phixpast1[xedge+ne_airx+1]) + 2*phiypast1[i] - phiypast2[i];
            elseif (colnum == ne_airx+1) # if the edge is at the right vertical interface
                phiy[i] = A2*(ly/lx)*(dt^2)/eps*ecurrent[real_ind] - (psiynow[real_ind-Ne_x]-psiypast1[real_ind-Ne_x])*dt - A1*(dt^2)/(mu*eps*lx^2)*(-phiypast1[i-1] + phiypast1[i] + phixpast1[xedge] - phixpast1[xedge+ne_airx]) + 2*phiypast1[i] - phiypast2[i];
            end
              
        end
    end
    return phix, phiy;
    
end