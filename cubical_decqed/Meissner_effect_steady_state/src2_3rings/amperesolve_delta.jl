function ampere_timedomain_delta(MMat::Array{Float64,2},Kcol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1})
    phip_past1 = [phixp_past1; phiyp_past1];
    phip_now   =  phip_past1 - inv(MMat)*Kcol;
    return phip_now;
end