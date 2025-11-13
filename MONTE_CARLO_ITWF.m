clear
clc
rng('default')
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
Realiz_num = 200;
N_f_vec = [8];
ric_fact = 0; 
MOD_CON_VEC = [0 1 2];
N_vec = [4];
Res_vec = zeros(length(MOD_CON_VEC), length(N_vec), length(N_f_vec), length(alph_chan_vec), 1);
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
            B_struct = ones(N, N);
            elem_num = N*(N + 1)/2;
            Z_0 = 50;
            Z = 1i.*ones(N);
            Z = (Z + Z.') ./ 2;
            
            %% MODE and Strcuture Selection
            if MODE_CON == 1
                for i = 1:N
                    for j = 1:N
                        if abs(i - j) > 1
                            B_struct(i, j) = 0;
                        end
                    end
                end
            end 
            if MODE_CON == 2
                for i = 1:N
                    for j = 1:N
                        if i~= j
                            B_struct(i, j) = 0;
                        end
                    end
                end
            end
            Z(B_struct == 0) = 0;
            [B_vec, ~] = vectorize_symmetric_matrix(B_struct);
            for nn = 1:length(N_f_vec)
                N_f = N_f_vec(nn);
                for realiz = 1:Realiz_num
                    [H_I, H_R, h_dir] = Rayleigh_channel_time(di, dr, N,N_f, L, path_exp, ref_loss_db, f0, BW, alph_chan);
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
                    delta_omeg = delta_ratio*(1/(N*norm(inv(Z + Z_0.*eye(N)), inf)));
                    omegaa = delta_omeg.*ones(N, N);
                    [w_vec, P] = vectorize_symmetric_matrix(omegaa);
                    cons_vec = zeros(N_f, 1);
                    P_n_vec = [];
                    Al = inv(Z_iter + Z_0.*eye(N));
                    Bl = Z_iter - Z_0.*eye(N);
                    for n = 1:N_f
                        anl = H_R(:, n).'*Al;
                        bnl = Bl*H_I(:, n);
                        cons_vec(n) = trace(anl*bnl);
                        Fn = Al*bnl*anl;
                        fn = reshape(Fn, [1 size(Fn, 1)*size(Fn, 2)]);
                        P_n_vec(n, :) = fn*P;
                    end
                    for n = 1:N_f
                        equiv_chan(n) = cons_vec(n) - P_n_vec(n, :)*w_vec;
                    end
                    h_var_0 = equiv_chan;
                    [coeff, f0value] = aprox_fixed_w(h_var_0, w0, K2, K4, N_f, 1, MODE_P);
                    f0 = -10;
                    fstar = 10;
                    f0_final = -10;
                    fstar_final = 10;
                    result_vec_alter = [];
                    result_vec_weights = [];
                    result_vec_phases = [];
                    result_vec_alter(end + 1) = real(f0value);
                    inc = 1;
                    out_iter = 1;
                    [coeff, ~] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, MODE_P);
                    lambdn_w = sqrt((sum(abs(coeff).^2))/(2*P_tx));
                    stall_counter = 0;
                    while norm(1 - fstar_final/f0_final) >= alt_thr
                        f0 = -10;
                        fstar = 10;
                        phase_iter = 1;
                        stall_counter_phase = 0;
                        lambdj_vec = ones(1, elem_num);
                        while norm(1 - fstar/f0) >= phase_thr
                            fstar = f0;
                            [coeff, ~] = aprox_fixed_w_dir(w_vec, w0, cons_vec, P_n_vec, K2, K4, N_f, MODE_P);
                            delta_omeg = inc*(delta_ratio)*(1/(norm(inv(Z_iter + Z_0.*eye(N)), inf)));
                            lambdj_vec = abs(coeff) / (2 * delta_omeg);
                            wstar = ((coeff)) ./ (2 * lambdj_vec);
                            w_vec = w_vec + up_delt_Z.*(wstar.' - w_vec);
                            w_vec(B_vec == 0) = 0;
                            omegaa = reshape(P*w_vec, [N N]);
                            Z_iter = Z_iter + 1i.*imag(omegaa);
                            w_vec = vectorize_symmetric_matrix_w(omegaa);
                            cons_vec = zeros(N_f, 1);
                            P_n_vec = [];
                            Al = inv(Z_iter + Z_0.*eye(N));
                            Bl = Z_iter - Z_0.*eye(N);
                            for n = 1:N_f
                                anl = H_R(:, n).'*Al;
                                bnl = Bl*H_I(:, n);
                                cons_vec(n) = trace(anl*bnl);
                                Fn = Al*bnl*anl;
                                fn = reshape(Fn, [1 size(Fn, 1)*size(Fn, 2)]);
                                P_n_vec(n, :) = fn*P;
                            end
                            for n = 1:N_f
                                equiv_chan(n) = cons_vec(n) - P_n_vec(n, :)*w_vec;
                            end
                            h_var_0 = equiv_chan;
                            [coeff, ~] = aprox_fixed_w_dir(w_vec, w0, cons_vec, P_n_vec, K2, K4, N_f, MODE_P); 
                            f0 = real((coeff*w_vec));
                            result_vec_phases(out_iter, phase_iter) = real(f0);
                            if real(f0) < max(result_vec_phases(out_iter, :))
                                stall_counter_phase = stall_counter_phase + 1;
                            end
                            if stall_counter_phase == stall_max_count_phase
                                break
                            end
                            phase_iter = phase_iter + 1;
                            if phase_iter == iter_max_in
                                break
                            end
                        end
                        PS = inv(Z_iter + Z_0.*eye(N))*(Z_iter - Z_0.*eye(N));
                        equiv_chan = zeros(N_f, 1);
                        for n = 1:N_f
                            equiv_chan(n) = H_R(:, n).'*PS*H_I(:, n);
                        end
                        equiv_chan = h_var_0;
                        f0 = -10;
                        fstar = 10;
                        weight_iter = 1;
                        stall_counter_weight = 0;
                        bet = 1;
                        w0 = zeros(N_f, 1);
                        for n = 1:N_f
                            w0(n) = exp(-1i*angle(equiv_chan(n)))*...
                                (norm(equiv_chan(n))^bet)*(sqrt((2*P_tx)/(sum(abs(equiv_chan).^(2*bet)))));
                        end
                        while norm(1 - fstar/f0) >= w_thr
                            fstar = f0;
                            [coeff, ~] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, MODE_P);
                            lambdn_w = sqrt((sum(abs(coeff).^2))/(2*P_tx));
                            wstar = coeff./lambdn_w;
                            absw0 = abs(w0) + up_delt_W.*(wstar - abs(w0));
                            w0 = absw0.*exp(-1i.*angle(equiv_chan));
                            [~, f0value] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, MODE_P);
                            f0 = real(real(sum(coeff.*absw0)));
                            result_vec_weights(out_iter, weight_iter) = real(f0);
                            weight_iter = weight_iter + 1;
                            if weight_iter == iter_max_in 
                                break
                            end
                            if real(f0) < max(result_vec_weights(out_iter, :))
                                stall_counter_weight = stall_counter_weight + 1;
                            end
                            if stall_counter_weight == stall_max_count_weight
                                break
                            end
                        end
                        fstar_final = f0_final;
                        out_iter = out_iter + 1;
                        [~, f0value] = aprox_fixed_w(w0, equiv_chan, K2, K4, N_f, 1, MODE_P);
                        f0_final = real(f0value);
                        result_vec_alter(end + 1) = real(f0value);
                        inc = inc*UP_Decrement;
                        if out_iter == iter_max_out
                            break
                        end
                        if real(f0_final) < max(result_vec_alter)
                            stall_counter = stall_counter + 1;
                        else
                            stall_counter = 0;
                        end
                        if stall_counter == stall_max_count_out
                            break
                        end
                        if out_iter == iter_max_out
                            break
                        end
                    end
                    
                    PS = inv(Z_iter + Z_0.*eye(N))*(Z_iter - Z_0.*eye(N));
                    equiv_chan = zeros(N_f, 1);
                    for n = 1:N_f
                        equiv_chan(n) = H_R(:, n).'*PS*H_I(:, n);
                    end
                    w0_cell{nmode, nant, nn, alph_c, realiz} = w0;
                    ris_cell{nmode, nant, nn, alph_c, realiz} = Z_iter;
                    h_var_0 = equiv_chan;
                    [~, f0value] = aprox_fixed_w(w0, h_var_0, K2, K4, N_f, 1, MODE_P);
                    FINAL_VALUE = real(f0value);
                    Res_vec(nmode, nant, nn) = Res_vec(nmode, nant, nn, alph_c) + max(result_vec_alter)/Realiz_num;
                    Res_vec_tot(nmode, nant, nn, alph_c, realiz) = FINAL_VALUE;
                    disp(strcat('IT_Final_RIS_', string(N), '_MODE_Power', string(MODE_P), '_MODERIS_', string(MODE_CON), ...
                        '_sub_carrier_', string(N_f), '_realization_', string(realiz)))
                    filename = strcat('IT_Final_RIS_', string(N), '_MODE_Power', string(MODE_P), '_MODERIS_', string(MODE_CON), ...
                        '_sub_carrier_', string(N_f), '_fade_factor', string(alph_chan),'.mat');
                    save(filename);
                end
            end
        end
    end
end