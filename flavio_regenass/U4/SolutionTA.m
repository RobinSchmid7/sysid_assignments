clear; close all; clc;

Ts = 1;

% TF
z = tf('z', -1);
G = 0.1 * z/(z^4-2.2*z^3+2.42*z^2-1.87*z+0.7225)

H = 0.5*tf([1 -0.9], [1 -0.25], Ts)

%% Sim

N = 2^10;
e = randn(N, 1) * sqrt(0.01); % Gaussian

% using Gaussian Distr Noise
u = randn(N, 1);
u = u/max(abs(u)); % Scale to satisfy the constraint

% Random Binary Input signal
% u = idinput(N);
% Periodic Input signal
% u2 = repmat(idinput(N/4), 4, 1);

t_vec = Ts * (0:N-1);

v = lsim(H, e, t_vec); % Noise
y = lsim(G, u, t_vec) + v; % Noisy Output

figure;
plot(t_vec, u); hold on; plot(t_vec, y);
xlim([0, N-1]); legend({'u', 'y'}); grid on;

% b) Unsmoothed ETFE
Y_N = fft(y);
U_N = fft(u);
G_N = Y_N(1:N/2+1) ./ U_N(1:N/2+1); % Only pos. Freq

% Continuous time frequencies
omega = 2*pi/(Ts*N)*(0:N/2); % Half of feasible freq Interval

% True System Freq Resp
G_freqresp = squeeze(freqresp(G, omega));
H_freqresp = squeeze(freqresp(H, omega));

% Plot ETFE, true system and magnitude of errorrs
fig1 = figure(1);
subplot(2,1, 1);
set(fig1, 'defaulttextinterpreter', 'latex');
% etc...


% c)
N_av = N/4;
u_reshaped = reshape(u, N_av, 4);
y_reshaped = reshape(y, N_av, 4);

U_N_reshaped = fft(u_reshaped, [], 1);
Y_N_reshaped = fft(y_reshaped, [], 1);
G_N_split = Y_N_reshaped(1:N_av/2+1, :) ./ U_N_reshaped(1:N_av/2+1, :);
    
G_N_averaged = mean(G_N_split, 2);


