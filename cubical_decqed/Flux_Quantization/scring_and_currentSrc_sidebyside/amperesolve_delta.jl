using LinearAlgebra
using IterativeSolvers

function ampere_timedomain_delta3D(MMat::Array{Float64,2},Kcol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, phizp_past1::Array{Float64,1})
    phip_past1 = [phixp_past1; phiyp_past1; phizp_past1];
    BLAS.set_num_threads(4);
    phip_now   =  phip_past1 - (MMat \ Kcol);
    BLAS.set_num_threads(1);
    return phip_now;
end


function ampere_timedomain_delta3D_gmres(MMat::Array{Float64,2},Kcol::Array{Float64,2},phixp_past1::Array{Float64,1}, phiyp_past1::Array{Float64,1}, phizp_past1::Array{Float64,1}, reltol::Float64)
    phip_past1 = [phixp_past1; phiyp_past1; phizp_past1];
    phip_delta, history = gmres(MMat, Kcol, log=true, restart=20, maxiter=1000, reltol=reltol);
    phip_now   =  phip_past1 - phip_delta;
    return phip_now;
end
