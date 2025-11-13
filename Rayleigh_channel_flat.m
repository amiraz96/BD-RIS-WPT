function [H_I, H_R, h_dir] = Rayleigh_channel_flat(di, dr, N, N_f, L, path_exp, ref_loss_db)


ref_loss = db2pow(ref_loss_db);
PL_i = di^(-path_exp)*ref_loss;
PL_r = dr^(-path_exp)*ref_loss;
tetang = pi/1.1;
d_dir = sqrt(di^2 + dr^2 - 2*di*dr*cos(tetang));


H_I_temp = sqrt(PL_i).*(1/sqrt(2)).*(randn(N, 1) + 1i*randn(N, 1));
H_R_temp = sqrt(PL_r).*(1/sqrt(2)).*(randn(N, 1) + 1i*randn(N, 1));

for n = 1:N_f
    H_I(:, n) = H_I_temp;
    H_R(:, n) = H_R_temp;
end

h_dir = zeros(N_f, 1);

end