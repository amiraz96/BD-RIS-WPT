function [H_I, H_R, h_dir] = Rayleigh_channel(di, dr, N, N_f, L, path_exp, ref_loss_db)


ref_loss = db2pow(ref_loss_db);
PLi = (di^(-path_exp))*ref_loss;
PLr = (dr^(-path_exp))*ref_loss;

rho = rand(1, L); % Random power distribution
rho = rho / sum(rho); % Normalize so that sum(rho) = 1
H_I_temp = zeros(N, N_f, L);
for l = 1:L
    H_I_temp(:,:,l) = sqrt(rho(l)) * (randn(N, N_f) + 1i * randn(N, N_f)) / sqrt(2);
end
H_I = sum(H_I_temp, 3);
H_I = sqrt(PLi).*H_I;


rho = rand(1, L); % Random power distribution
rho = rho / sum(rho); % Normalize so that sum(rho) = 1
H_R_temp = zeros(N, N_f, L);
for l = 1:L
    H_R_temp(:,:,l) = sqrt(rho(l)) * (randn(N, N_f) + 1i * randn(N, N_f)) / sqrt(2);
end
H_R = sum(H_R_temp, 3);
H_R = sqrt(PLr).*H_R;
h_dir = zeros(N_f, 1);

end