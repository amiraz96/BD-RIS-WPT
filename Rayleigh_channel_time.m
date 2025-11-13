function [H_I, H_R, h_dir] = Rayleigh_channel_time(di, dr, N,N_f, L, path_exp, ref_loss_db, f0, BW, alph_chan)


ref_loss = db2pow(ref_loss_db);
PLi = (di^(-path_exp))*ref_loss;
PLr = (dr^(-path_exp))*ref_loss;
% Calculate the coherence bandwidth
Bc = alph_chan * BW;
% Calculate the RMS delay spread (tau_rms) based on the coherence bandwidth
tau_rms = 1 / Bc;
% Delay spread for the channel
tau_max = 2 * tau_rms; % Assume maximum delay spread is 2 times the RMS delay spread
% Subcarrier spacing
delta_f = BW / N_f;  % Subcarrier spacing (Hz)
% Generate the delay values (uniformly spaced up to tau_max)
delays = linspace(0, tau_max, L);
% Generate random Rayleigh fading for each tap and each antenna
% This models the multipath components with Rayleigh fading
% Generate random power values for each tap
tap_powers = rand(1, L);
% Normalize the tap powers so that they sum to 1
tap_powers = tap_powers / sum(tap_powers);
% Generate random Rayleigh fading for each tap and each antenna
% This models the multipath components with Rayleigh fading
h_taps = (randn(N, L) + 1i * randn(N, L)) / sqrt(2);
% Scale the taps by the square root of the normalized power values
h_taps = h_taps .* sqrt(tap_powers);
% Initialize channel frequency response matrix
H_I = zeros(N, N_f);
% Calculate the channel frequency response for each subcarrier
for k = 1:N_f
    % Frequency corresponding to subcarrier k, shifted by the carrier frequency f0
    f_k = f0 + (k - 1 - N_f/2) * delta_f;
    
    % Frequency response for each tap and sum them up
    H_I(:, k) = sum(h_taps .* exp(-1i * 2 * pi * f_k * delays), 2);
end
H_I = H_I * sqrt(PLi);


tap_powers = rand(1, L);
% Normalize the tap powers so that they sum to 1
tap_powers = tap_powers / sum(tap_powers);
% Generate random Rayleigh fading for each tap and each antenna
% This models the multipath components with Rayleigh fading
h_taps = (randn(N, L) + 1i * randn(N, L)) / sqrt(2);
% Scale the taps by the square root of the normalized power values
h_taps = h_taps .* sqrt(tap_powers);
% Initialize channel frequency response matrix
H_R = zeros(N, N_f);
% Calculate the channel frequency response for each subcarrier
for k = 1:N_f
    % Frequency corresponding to subcarrier k, shifted by the carrier frequency f0
    f_k = f0 + (k - 1 - N_f/2) * delta_f;
    
    % Frequency response for each tap and sum them up
    H_R(:, k) = sum(h_taps .* exp(-1i * 2 * pi * f_k * delays), 2);
end
H_R = H_R * sqrt(PLr);


h_dir = zeros(N_f, 1);


end