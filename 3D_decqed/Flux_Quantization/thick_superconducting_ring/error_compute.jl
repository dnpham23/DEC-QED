using Statistics

Bfield_error = zeros(Nstep*nperiod);
phix_error  = zeros(Nstep*nperiod);
phiy_error  = zeros(Nstep*nperiod);
psix_error   = zeros(Nstep*nperiod);
psiy_error   = zeros(Nstep*nperiod);
Bcount = 0;
phix_count = 0;
phiy_count = 0;
psix_count  = 0;
psiy_count  = 0;
temp = 0.0;
#for tn = 1:Nstep*nperiod
for tn = 1:Nstep*nperiod
    temp = 0.0;
    Bcount = 0;
    for i = 1:ne_x
        for j = 1:ne_y
            if (Bgrid[i,j,tn] >= 0.001)
                temp += abs((Bgrid[i,j,tn] - Bgrid1[i,j,tn])/Bgrid[i,j,tn]);
                Bcount += 1;
            end
        end
    end
    Bfield_error[tn] = temp/Bcount;
    
    temp = 0.0;
    phiy_count = 0;
    for i = 1:Nx
        for j = 1:ne_y
            if (phiygrid[i,j,tn] != 0)
                temp += abs((phiygrid[i,j,tn] - phiygrid1[i,j,tn])/phiygrid[i,j,tn]);
                phiy_count += 1;
            end
        end
    end
    phiy_error[tn] = temp/phiy_count;
end