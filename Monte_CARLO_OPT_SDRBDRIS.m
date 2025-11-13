clear
clc
% cvx_solver mosek;
rng(25)
f0 = 2.4e9;
lamb_f = 3e8/f0;
d = 2;
di = 2;
dr = 2;
ref_loss_db = -20*log10((4*pi*1)/lamb_f);
L = 18;
path_exp = 2;
K2 = 0.17;
K4 = 957.25;
P_tx_dB = 50; %in dBm
P_tx = (10^((P_tx_dB - 30)/10));
BW = 10e6;

%% OPT PARAMS
delta_ratio = 0.01;
phase_thr = 1e-6;
w_thr = 1e-6;
alt_thr = 1e-4;
UP_Decrement = 1;
up_delt_Z = 1;
up_delt_W = 0.6;
iter_max_out = 100;
iter_max_in = 500;
stall_max_count_phase = 20;
stall_max_count_weight = 10;
stall_max_count_out = 5;


alph_chan_vec = [0.1];
Realiz_num = 1;
N_f_vec = [8];
ric_fact = 0; 
MOD_CON_VEC = [0];
N_vec = [4];
Res_vec = zeros(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), 1);
Res_vec_obj = zeros(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), 1);
Res_vec_tot = zeros(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), Realiz_num);
MODE_P = 1; % 1 for DC and 2 for RF
w0_cell = cell(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), Realiz_num);
ris_cell = cell(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), Realiz_num);

for alph_c = 1:length(alph_chan_vec) 
    alph_chan = alph_chan_vec(alph_c);
    for nmode = 1:length(MOD_CON_VEC)
        MODE_CON = MOD_CON_VEC(nmode); %0 for fully, 1 for Tricon, 2 for Diag
        for nant = 1:length(N_vec)
            N = N_vec(nant);
            elem_num = N*(N + 1)/2;
            Z_0 = 50;
            Z = 1i.*ones(N);
            Z = (Z + Z.') ./ 2;
            for nn = 1:length(N_f_vec)
                N_f = N_f_vec(nn);
                w0 = rand(N_f, 1) + 1i.*rand(N_f, 1)./10000;
                v_vec = [];
                H = [];
                a_vecs = [];
                for realiz = 1:Realiz_num
                    [H_I, H_R, h_dir] = Rayleigh_channel_time(di, dr, N,N_f, L, path_exp, ref_loss_db, f0, BW, alph_chan);
%                     [H_I, H_R, h_dir, H_los, H_R1] = Rician_channel(di, dr, N, N_f, L, path_exp, ref_loss_db, ric_fact, f0);
                    PS = inv(Z + Z_0.*eye(N))*(Z - Z_0.*eye(N));
                    equiv_chan = zeros(N_f, 1);
                    for n = 1:N_f
                        equiv_chan(n) = H_R(:, n).'*PS*H_I(:, n);
                    end
                    bet = 3;
                    w0 = zeros(N_f, 1);
                    for n = 1:N_f
                        w0(n) = exp(-1i*angle(equiv_chan(n)))*...
                            (norm(equiv_chan(n))^bet)*(sqrt((2*P_tx)/(sum(abs(equiv_chan).^(2*bet)))));
                    end
                    Z_iter = Z;
                    [phivec, P] = vectorize_symmetric_matrix(PS);
                    for n = 1:N_f
                        H(:, :, n) = H_I(:, n)*H_R(:, n).';
                    end
                    
                    for n = 1:N_f
                        a_vecs(n, :) = w0(n).*reshape(H(:, :, n), [1, N^2])*P;
                    end
                    D = zeros(N_f, elem_num, elem_num);
                    intval = 0;
                    for nout = 1:N_f
                        for n = 1:N_f
                            if n + intval > N_f
                                break
                            end
                            D(nout, :, :) = D(nout, :, :) + reshape(a_vecs(n, :)'*a_vecs(n + intval, :), [1 elem_num elem_num]);
                        end
                        intval = intval + 1;
                    end
                    kvec = ones(N_f, 1).*(3/4)*K4;
                    kvec(1) = (3/8)*K4;
                    K0 = diag(kvec);
                    dvec = zeros(N_f, 1);
                    for n = 1:N_f
                        Dk = reshape(D(n, :, :), [elem_num elem_num]);
                        X = phivec*phivec';
                        dvec(n) = trace(Dk*X);
                    end
                    fout = 1000;
                    for out_iter = 1:100
                        fin = 1000;
                        for scaiter = 1:100
                            D0 = reshape(D(1, :, :), [elem_num elem_num]);
                            myJ = -(K2/4)*D0 - (3*K4/8)*dvec(1)*D0;
                            for n = 2:N_f
                                Dk = reshape(D(n, :, :), [elem_num elem_num]);
                                X = phivec*phivec';
                                dvec(n) = trace(Dk*X);
                                myJ = myJ - (3*K4/4)*dvec(n)'*Dk;
                            end
                            K1 = myJ + myJ';
                            cvx_begin sdp quiet
                                variable X(elem_num, elem_num) hermitian semidefinite
                                minimize trace(K1*X)
                                subject to 
                                    for i = 1:N
                                        for j = i:N
                                            Pi = P((i-1)*N + (1:N), :);
                                            Pj = P((j-1)*N + (1:N), :);
                                            Pij = Pi'*Pj;
                                            if i == j
                                                trace(Pij*X) == 1;
                                            else
                                                trace(Pij*X) == 0;
                                            end
                                        end
                                    end
                            cvx_end
                            sca_res(out_iter, scaiter) = trace(K1*X)
                            [V, F] = eigs(X);
                            eigenvalues = diag(F);
                            [~, idx] = max(eigenvalues);
                            dominant_eigenvector = V(:, idx);
                            dominant_eigenvalue = eigenvalues(idx);
                            phivec = sqrt(dominant_eigenvalue) * dominant_eigenvector;
                            dvec = zeros(N_f, 1);
                            for n = 1:N_f
                                Dk = reshape(D(n, :, :), [elem_num elem_num]);
                                dvec(n) = trace(Dk*X);
                            end
                            myJ = -(K2/4)*D0 - (3*K4/8)*dvec(1)*D0;
                            for n = 2:N_f
                                Dk = reshape(D(n, :, :), [elem_num elem_num]);
                                dvec(n) = trace(Dk*X);
                                myJ = myJ - (3*K4/4)*dvec(n)'*Dk;
                            end
                            if norm(1 - trace(K1*X)/fin) <= 1e-5
                                break
                            end
                            fin = trace(K1*X);
                            K1 = myJ + myJ';
                        end
                        phivec = gaussian_randomization_cholesky(X, 50000, K1);
                        PS = reshape(P*phivec, [N N]);
                        equiv_chan = zeros(N_f, 1);
                        for n = 1:N_f
                            equiv_chan(n) = H_R(:, n).'*PS*H_I(:, n) + h_dir(n);
                        end
                        w0 = zeros(N_f, 1);
                        bet = 3;
                        f0 = -10;
                        fstar = 10;
                        weight_iter = 1;
                        w_thr = 1e-6;
                        stall_counter_weight = 0;
                        up_delt_W = 0.5;
                        iter_max_in = 100;
                        stall_max_count = 10;
                        for n = 1:N_f
                            w0(n) = exp(-1i*angle(equiv_chan(n)))*...
                                (norm(equiv_chan(n))^bet)*(sqrt((2*P_tx)/(sum(abs(equiv_chan).^(2*bet)))));
                        end
                        while norm(1 - fstar/f0) >= w_thr
                            fstar = f0;
                            [coeff, ~] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, 1);
                            lambdn_w = sqrt((sum(abs(coeff).^2))/(2*P_tx));
                            wstar = coeff./lambdn_w;
                            absw0 = abs(w0) + up_delt_W.*(wstar - abs(w0));
                            w0 = absw0.*exp(-1i.*angle(equiv_chan));
                            [~, f0value] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, 1);
                            f0 = real(real(sum(coeff.*absw0)));
                            weight_iter = weight_iter + 1;
                            if weight_iter == iter_max_in 
                                break
                            end
                        end
                        [~, f0value] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, 1);
                        for n = 1:N_f
                            a_vecs(n, :) = w0(n).*reshape(H(:, :, n), [1, N^2])*P;
                        end
                        D = zeros(N_f, elem_num, elem_num);
                        intval = 0;
                        for nout = 1:N_f
                            for n = 1:N_f
                                if n + intval > N_f
                                    break
                                end
                                D(nout, :, :) = D(nout, :, :) + reshape(a_vecs(n, :)'*a_vecs(n + intval, :), [1 elem_num elem_num]);
                            end
                            intval = intval + 1;
                        end
                        D0 = reshape(D(1, :, :), [elem_num elem_num]);
                        dvec = zeros(N_f, 1);
                        for n = 1:N_f
                            Dk = reshape(D(n, :, :), [elem_num elem_num]);
                            dvec(n) = trace(Dk*X);
                        end
                        final_val = 0.5*K2*phivec'*D0*phivec + dvec'*K0*dvec;
                        newval = 0.5*K2*phivec'*D0*phivec + (3/8)*K4*phivec'*D0*phivec*(phivec'*D0*phivec)';
                        for n = 2:N_f
                            Dk = reshape(D(n, :, :), [elem_num elem_num]);
                            newval = newval + (3/4)*K4*phivec'*Dk*phivec*(phivec'*Dk*phivec)';
                        end
                        optval = 0.5*K2*trace(X*D0) + (3/8)*K4*trace(X*D0*X*D0');
                        for n = 2:N_f
                            Dk = reshape(D(n, :, :), [elem_num elem_num]);
                            optval = optval + (3/4)*K4*trace(X*Dk*X*Dk');
                        end
                        res_alter(out_iter) = optval
                        if norm(1 - optval/fout) <= 1e-5
                            break
                        end
                        fout = optval;
                    end
                    [U, ~] = svd(PS);
                    PS = takagi_phase_randomization_approximation(PS, 500000, H_R, H_I, w0, K2, K4);
                    equiv_chan = zeros(N_f, 1);
                    for n = 1:N_f
                        equiv_chan(n) = H_R(:, n).'*PS*H_I(:, n) + h_dir(n);
                    end
                    [~, f0value] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, 1);
                    FINAL_VALUE = real(f0value);
                    FINAL_VALUE_OPT = real(optval);
                    Res_vec(nmode, nant, nn, alph_c) = Res_vec(nmode, nant, nn, alph_c) + real(f0value)/Realiz_num;
                    Res_vec_obj(nmode, nant, nn, alph_c) = Res_vec_obj(nmode, nant, nn, alph_c) + real(optval)/Realiz_num;
                    Res_vec_tot(nmode, nant, nn, alph_c, realiz) = real(f0value);
                    disp(strcat('SDR_Final_RIS_', string(N), '_MODE_Power', string(MODE_P), '_MODERIS_', string(MODE_CON), ...
                        '_sub_carrier_', string(N_f), '_realization_', string(realiz)))
                    filename = strcat('SDR_Final_RIS_', string(N), '_MODE_Power', string(MODE_P), '_MODERIS_', string(MODE_CON), ...
                        '_sub_carrier_', string(N_f), '_fade_factor', string(alph_chan),'_ricfact_', string(ric_fact),'.mat');
%                     save(filename);
                end
            end
        end
    end
end




