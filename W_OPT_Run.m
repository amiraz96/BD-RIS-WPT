function [w0, f0] = W_OPT_Run(w0, h0, N_f, K2, K4, obj_scale)
P_tx  = 1;

cvx_begin quiet
    variable w(N_f, 1) complex
    [coeff, f0value] = aprox_fixed_w(w0, h0, K2, K4, N_f, 1);
    objval = 0;
    for nin = 1:N_f
        objval = objval + (coeff(nin)*(w(nin)));
    end
    maximize(obj_scale.*real(objval))
    subject to
        0.5.*w'*w <= P_tx
cvx_end
w0 = w;
f0 = real(objval);

end