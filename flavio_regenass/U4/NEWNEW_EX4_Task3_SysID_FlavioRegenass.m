%-------------------------------------------------
% System Identification Exercise Set 4, Task 3
% Flavio Regenass
% October 2020
%-------------------------------------------------
% Empirical Transfer Function Estimate ETFE

clear; close all; clc;
% Sampling Time
Ts = 1;

% TF
z = tf('z', -1);
G = 0.1 * z/(z^4-2.2*z^3+2.42*z^2-1.87*z+0.7225)

H = 0.5*tf([1 -0.9], [1 -0.25], Ts)

%% Sim

N = 2^10;
e = randn(N, 1) * sqrt(0.01); % Gaussian

% using Gaussian Distr Noise
% u = randn(N, 1);
% u = u/max(abs(u)); % Scale to satisfy the constraint

% Random Binary Input signal
u = idinput(N);
% Periodic Input signal
u2 = repmat(idinput(N/4), 4, 1);

% Timesteps
t_vec = Ts * (0:N-1);

% Simulation
v = lsim(H, e, t_vec); % Noise
y = lsim(G, u, t_vec) + v; % Noisy Output

figure;
plot(t_vec, u); hold on; plot(t_vec, y);
xlim([0, N-1]); legend({'u', 'y'}); grid on;

% Simulation for Periodic Input
v2 = lsim(H, e, t_vec); % Noise
y2 = lsim(G, u2, t_vec) + v2; % Noisy output

figure;
plot(t_vec, u2); hold on; plot(t_vec, y2);
xlim([0, N-1]); legend({'u', 'y'}); grid on;

%% b) Unsmoothed ETFE
Y_N = fft(y); 
U_N = fft(u); % We cannot do this with Periodic Signals, since we would receive value 0 for all disharmonic frequencies ...
% which would result in division by zero. So we must Do the averaged procedure!
G_N = Y_N(1:N/2+1) ./ U_N(1:N/2+1); % Only pos. Frequencies. Our FFT is symmetrical since the input signal is real valued.

% Continuous time frequencies to sample System Freq Response
omega = 2*pi/(Ts*N)*(0:N/2); % Half of feasible freq Interval

% True System Freq Resp
G_freqresp = squeeze(freqresp(G, omega));
H_freqresp = squeeze(freqresp(H, omega));

% Plot ETFE, true system and magnitude of errorrs
fig1 = figure(1);
subplot(2,1, 1);
set(fig1, 'defaulttextinterpreter', 'latex');
loglog(omega, abs(G_N), 'b--', omega, abs(G_freqresp), 'g', ...
    omega, abs(H_freqresp), 'r');
xlim([0.2, max(omega)]);
ylim([1e-2, 1e3]);
title('Unsmoothed ETFE');
xlabel('Frequency rad(sample)');
ylabel('Magnitude');
legend('ETFE unsmoothed: |G_N(z)|', 'True Plant: |G(z)|', ...
    'Noise System |H(z)|');

subplot(2, 1, 2)
loglog(omega, abs(G_N - G_freqresp), 'b--', omega, abs(H_freqresp), 'r');
xlim([0.2, max(omega)]);
ylim([1e-2, 1e3]);
title('Magnitude of Errors');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
legend('Error: |G_N(z) - G(z)|', 'Noise System |H(z)|');


%% c) Averaged ETFE
%Split Data into 4 Parts
N_av = N/4;
u_reshaped = reshape(u, N_av, 4); % Split Data into 4 Columns
y_reshaped = reshape(y, N_av, 4);
u_reshaped2 = reshape(u2, N_av, 4);
y_reshaped2 = reshape(y2, N_av, 4);

% FFT on each record separately
U_N_reshaped = fft(u_reshaped, [], 1); % Specifiy Dimension, here we want fft for each column! (Dimension 1)
Y_N_reshaped = fft(y_reshaped, [], 1);
U_N_reshaped2 = fft(u_reshaped2, [], 1);
Y_N_reshaped2 = fft(y_reshaped2, [], 1);

% ETFE of each record
G_N_split = Y_N_reshaped(1:N_av/2+1, :) ./ U_N_reshaped(1:N_av/2+1, :);
G_N_split2 = Y_N_reshaped2(1:N_av/2+1, 2:4) ./ U_N_reshaped2(1:N_av/2+1, 2:4);
% Neglect the first set, since this contains the transient response,
% which we do not want to include. We can do this only since we measure a
% periodic signal!

% Average over all 4 ETFE's (along dimension 2)
G_N_averaged = mean(G_N_split, 2); % Mean of all four G_N_split at each calculated frequency
G_N_averaged2 = mean(G_N_split2, 2);

omega_av = 2*pi/(Ts * N_av) * (0:N_av/2); % Through using averaged ETFE, our frequency resolution becomes smaller

% Plot Averaged ETFE
fig2 = figure(3);
set(fig2, 'defaulttextinterpreter', 'latex');
loglog(omega_av, abs(G_N_averaged), 'b--', omega, abs(G_N), 'g', ...
    omega, abs(G_freqresp), 'r');
hold on; loglog(omega_av, abs(G_N_averaged2), 'k--');
xlim([0.2, max(omega_av)]);
ylim([1e-2, 1e3]);
title('Averaged ETFE');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
% legend('ETFE averaged: |G_{av}(z)|', 'ETFE unsmoothed: |G_N(z)|', 'True plant: |G(z)|');
legend('ETFE averaged: |G_{av}(z)|', 'ETFE unsmoothed: |G_N(z)|', 'True plant: |G(z)|', 'Periodic Average: |G_{av2}(z)|');

