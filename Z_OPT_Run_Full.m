function [f0, h_var, omegaa] = Z_OPT_Run_Full(Z_iter, w0, h_var_0, H_I, H_R, ...
    N_f, N, K2, K4, delta_omeg, Z_0)

cvx_begin quiet
    variable h_var(N_f, 1) complex
    variable omegaa(N,N) complex symmetric
    [coeff, f0value] = aprox_fixed_w(h_var_0, w0, K2, K4, N_f, 1);
    objval = 0;
    for nin = 1:N_f
%         objval = objval + (coeff(nin)*(h_var(nin) - h_var_0(nin)));
        objval = objval + coeff(nin)*(h_var(nin));
    end
    maximize(real(objval))
    subject to
        Al = inv(Z_iter + Z_0.*eye(N));
        Bl = Z_iter - Z_0.*eye(N);
        for n = 1:N_f
            anl = H_R(:, n).'*Al;
            bnl = Bl*H_I(:, n);
            h_var(n, 1) == anl*(eye(N) - omegaa*Al)*bnl;
        end
        for nn = 1:N
            for nn2 = nn:N
                norm(omegaa(nn, nn2)) <= delta_omeg;
            end
        end
cvx_end

f0 = real(objval);


end