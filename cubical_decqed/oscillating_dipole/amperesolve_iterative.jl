using LinearAlgebra

function ampere_timedomain_newton(AMat::Array{Float64,2},BMat::Array{Float64,2},Ccol::Array{Float64,2},iguess::Array{Float64,2}, Ne_sc::Int64,max_iter::Int64,conv_tol::Float64)

    phip_prev = iguess;
    phip_next = zeros(length(iguess),1);
    phip_diag = Diagonal(vec(phip_prev));
    mean_err = 0;
    err = zeros(length(iguess),1);
    for i = 1:max_iter
        J = 2*AMat*phip_diag + BMat;
        phip_next = phip_prev - inv(J)*(AMat*(phip_prev.*phip_prev) + BMat*phip_prev + Ccol);
        err = abs.((phip_next-phip_prev)./phip_prev);
        mean_err = sum(err)/length(err);
        if mean_err <= conv_tol
             break;
        end
        phip_prev = phip_next;
        phip_diag = Diagonal(vec(phip_prev));
    end
    if mean_err > conv_tol
       print("solver didn't converge after ",max_iter," iterations\n");
    end
    return phip_next, mean_err, err;
end


