function [H_I, H_R, h_dir, los_comp, H_R1] = Rician_channel(di, dr, N, N_f, L, path_exp, ref_loss_db, ric_fact, f_0)

ref_loss = db2pow(ref_loss_db);
PLi = (di^(-path_exp))*ref_loss;
PLr = (dr^(-path_exp))*ref_loss;
tetang = pi/1.1;
Nvecmain = [4 8 12 16 20 24 28 32 40 64 128];
Nvec = [2 2; 2 4; 3 4; 4 4; 4 5; 4 6; 4 7; 4 8; 5 8; 8 8; 8 16];

for nn = 1:length(Nvecmain)
    if N == Nvecmain(nn)
        Ni = Nvec(nn, 1);
        Nh = Nvec(nn, 2);
    end
end
LOS_fact = sqrt(ric_fact/ (ric_fact + 1));
NLOS_fact = sqrt(1/ (ric_fact + 1));
theta_i = pi/8;            % Azimuth angle of arrival (in radians)
phi_i = pi/8;              % Elevation angle of arrival (in radians)
lambdaa = (3e8/f_0);
d_int = lambdaa/2;
theta_r = pi/8;            % Azimuth angle of arrival (in radians)
phi_r = pi/8;              % Elevation angle of arrival (in radians)

rho = rand(1, L); % Random power distribution
rho = rho / sum(rho); % Normalize so that sum(rho) = 1
H_I_temp = zeros(N, N_f, L);
for l = 1:L
    H_I_temp(:,:,l) = sqrt(rho(l)) * (randn(N, N_f) + 1i * randn(N, N_f)) / sqrt(2);
end
los_comp = zeros(N, N_f);
for m = 1:Ni
    for n = 1:Nh
        % Phase shift for the (m,n)th element
        delta_psi_mn = (2 * pi / lambdaa) * ...
                       ((m - 1) * d_int * sin(theta_i) * cos(phi_i) + ...
                        (n - 1) * d_int * sin(theta_i) * sin(phi_i));
        % LOS component for the (m,n)th element
        los_comp((m - 1)*Nh + n, :) = ones(1, N_f).*exp(1i * delta_psi_mn);
    end
end
H_I = sum(H_I_temp, 3);
H_I = sqrt(PLi).*(NLOS_fact.*H_I + LOS_fact.*los_comp);


rho = rand(1, L); % Random power distribution
rho = rho / sum(rho); % Normalize so that sum(rho) = 1
H_R_temp = zeros(N, N_f, L);
for l = 1:L
    H_R_temp(:,:,l) = sqrt(rho(l)) * (randn(N, N_f) + 1i * randn(N, N_f)) / sqrt(2);
end
los_comp = zeros(N, N_f);
for m = 1:Ni
    for n = 1:Nh
        % Phase shift for the (m,n)th element
        delta_psi_mn = (2 * pi / lambdaa) * ...
                       ((m - 1) * d_int * cos(theta_r) * cos(phi_r) + ...
                        (n - 1) * d_int * cos(theta_r) * sin(phi_r));
        % LOS component for the (m,n)th element
        los_comp((m - 1)*Nh + n, :) = ones(1, N_f).*exp(1i * delta_psi_mn);
    end
end
H_R1 = sum(H_R_temp, 3);
H_R = sqrt(PLr).*(NLOS_fact.*H_R1 + LOS_fact.*los_comp);


% d_dir = sqrt(di^2 + dr^2 - 2*di*dr*cos(tetang));
% PL_dir = (d_dir^(-path_exp))*ref_loss;
% rho = rand(1, L); % Random power distribution
% rho = rho / sum(rho); % Normalize so that sum(rho) = 1
% h_dir = zeros(N_f, L);
% for l = 1:L
%     h_dir(:,l) = sqrt(rho(l)) * (randn(N_f, 1) + 1i * randn(N_f, 1)) / sqrt(2);
% end
% los_comp = ones(N_f, 1);
% h_dir = sum(h_dir, 2);
% h_dir = PL_dir.*(NLOS_fact.*h_dir + LOS_fact.*los_comp);
h_dir = zeros(N_f, 1);
end