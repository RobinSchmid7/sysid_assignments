% System Identification Ex 10
% Robin Schmid, schmirob@ethz.ch
%% 1.
clear all; close all; clc;

% Noise variance
lambda = 0.1;

% Continuous time transfer function
zeta_z = 0.1;
omega_z = 3.0;
zeta_p = 0.1;
omega_p = 3.5;
s = tf('s');
G_s = (s^2 + 2*zeta_z*omega_z*s + omega_z^2)/(s^2 + 2*zeta_p*omega_p*s + omega_p^2) * 5000/((s+50)*(s+200));
%pole(G_s) % All poles have negative real part, continuous time O.L. is stable

% Discretize with zero order hold
Ts = 0.02;
G_d = c2d(G_s, Ts, 'zoh');
%abs(pole(G_d)); % All poles have magnitude smaller than one, discrete time O.L. is stable after sampling

% Controller
z = tf('z', Ts);
C_d = (1.25*z - 0.75)/(z - 1);

% Check stability of closed loop system
% Option 1: Check with Nyquist theorem, here no encirclements around -1 so
% the system is closed loop stable
% nyquist(G_d*C_d);
% Option 2: Check poles of closed loop system from closed loop transfer
% function T or equivalently S (have same poles)
S_d = 1/(1 + G_d*C_d);
abs(pole(S_d)); % All poles have magnitude smaller than one, discrete time C.L. is stable

%% 2.
% ETFE for G = Y/U
% PRBS reference and noise
M = 5; % Use 10 periods
N = 10;
K = 2^N-1; % Use 1023 data points per period, 2^x-1 s.t. prbs is periodic

r_period = 0.1*idinput(K);
r = repmat(r_period, M, 1);

v = sqrt(lambda)*randn(M*K,1);

t = Ts*(0:(M*K-1));

y = lsim(S_d*G_d*C_d, r, t) + lsim(S_d, v, t);
u = lsim(S_d*C_d, r, t) - lsim(S_d*C_d, v, t);

% Discard 1st period to improve transient error
y(1:K) = [];
u(1:K) = [];

% Time domain averaging to reduce effect of noise, filter over 4 periods
y_reshaped = reshape(y, K, 4);
u_reshaped = reshape(u, K, 4);

% Average over remaining periods in time domain
y_avg = mean(y_reshaped, 2);
u_avg = mean(u_reshaped, 2);

% ETFE of averaged input and output
Y_avg = fft(y_avg);
U_avg = fft(u_avg);

G_ETFE = Y_avg./U_avg;

% Frequencies, consider only positive ones
omega = (2*pi/(K*Ts))*(0:K-1);
idx = find(omega > 0 & omega < pi/Ts);

% True frequency response
G_freq = squeeze(freqresp(G_d, omega));

% Plotting
f1 = figure(1);
set(f1, 'visible', 'on');
% Magnitude
subplot(2,1,1)
loglog(omega(idx), abs(G_ETFE(idx)), 'r-');
hold on
loglog(omega(idx), abs(G_freq(idx)), 'b-');
title({'Magnitude of ETFE of Y and U and true frequency response'});
legend({'ETFE of', 'true frequency response G'}, 'Location', 'northwest');
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on
% Error magnitude
subplot(2,1,2)
loglog(omega(idx), abs(G_ETFE(idx)-G_freq(idx)), 'k-');
title({'Error of ETFE of Y and U'});
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on

% For high frequencies quite noisy

%% 3.
% Again use 5 periods and discard the first period
r_period = 0.1*idinput(K);
r = repmat(r_period, M, 1);

t = Ts*(0:(M*K-1));

% Measuring error
e = lsim(S_d, r, t) - lsim(S_d, v, t);

% ETFE for sensitivity S = E/R
% Discard 1st period to reduce transient error
e(1:K) = [];
r(1:K) = [];

% Time domain averaging to reduce effect of noise, filter over 4 periods
e_reshaped = reshape(e, K, 4);
r_reshaped = reshape(r, K, 4);

% Average over remaining periods in time domain
e_avg = mean(e_reshaped, 2);
r_avg = mean(r_reshaped, 2);

% ETFE of averaged input and output
E_avg = fft(e_avg);
R_avg = fft(r_avg);

S_ETFE = E_avg./R_avg;

% Frequencies, consider only positive ones
omega = (2*pi/(K*Ts))*(0:K-1);
idx = find(omega > 0 & omega < pi/Ts);

% True frequency response of sensitivity
S_freq = squeeze(freqresp(S_d, omega));

% Plotting
f2 = figure(2);
set(f2, 'visible', 'on');
% Magnitude
subplot(2,1,1)
loglog(omega(idx), abs(S_ETFE(idx)), 'r-');
hold on
loglog(omega(idx), abs(S_freq(idx)), 'b-');
title({'Magnitude of ETFE and true frequency response'});
legend({'ETFE of S', 'true frequency response of S'}, 'Location', 'northwest');
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on
% Error magnitude
subplot(2,1,2)
loglog(omega(idx), abs(S_ETFE(idx)-S_freq(idx)), 'k-');
title({'Error of ETFE of Y and U'});
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on

% Worse approximated at high frequencies

%% 4.
% Again use 5 periods and discard the first period
w_period = 0.1*idinput(K);
w = repmat(w_period, M, 1);

t = Ts*(0:(M*K-1));

% Measuring output
y = lsim(G_d*S_d, w, t) + lsim(S_d, v, t);

% ETFE for sensitivity G*S = Y/W
% Discard 1st period to reduce transient error
y(1:K) = [];
w(1:K) = [];

% Time domain averaging to reduce effect of noise, filter over 4 periods
y_reshaped = reshape(y, K, 4);
w_reshaped = reshape(w, K, 4);

% Average over remaining periods in time domain
y_avg = mean(y_reshaped, 2);
w_avg = mean(w_reshaped, 2);

% ETFE of averaged input and output
Y_avg = fft(y_avg);
W_avg = fft(w_avg);

GS_ETFE = Y_avg./W_avg;

% Frequencies, consider only positive ones
omega = (2*pi/(K*Ts))*(0:K-1);
idx = find(omega > 0 & omega < pi/Ts);

% True frequency response of G_d*S_d
GS_freq = squeeze(freqresp(G_d*S_d, omega));

% Plotting
f3 = figure(3);
set(f3, 'visible', 'on');
% Magnitude
subplot(2,1,1)
loglog(omega(idx), abs(GS_ETFE(idx)), 'r-');
hold on
loglog(omega(idx), abs(GS_freq(idx)), 'b-');
title({'Magnitude of ETFE and true frequency response'});
legend({'ETFE of G*S', 'true frequency response of G*S'}, 'Location', 'northwest');
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on
% Error magnitude
subplot(2,1,2)
loglog(omega(idx), abs(GS_ETFE(idx)-S_freq(idx)), 'k-');
title({'Error of ETFE of GS and S'});
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on

%% 5.
% Find G from estimates of G*S and S
G_combined = GS_ETFE./S_ETFE;

% Plotting
f4 = figure(4);
set(f4, 'visible', 'on');
% Magnitude
subplot(2,1,1)
loglog(omega(idx), abs(G_combined(idx)), 'r-');
hold on
loglog(omega(idx), abs(G_ETFE(idx)), 'b-');
loglog(omega(idx), abs(G_freq(idx)), 'k-');
title({'Magnitude ETFE of G*S and S, ETFE from Y and U and true frequency response'});
legend({'ETFE of G*S and S', 'ETFE of Y and U', 'true freq response'}, 'Location', 'northwest');
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on
% Error magnitude
subplot(2,1,2)
loglog(omega(idx), abs(G_combined(idx)-G_freq(idx)), 'r-');
hold on
loglog(omega(idx), abs(G_ETFE(idx)-G_freq(idx)), 'b-');
title({'Error of ETFE of G*S and S and ETFE of Y and U from true freq response'});
legend({'Error ETFE of G*S and S', 'Error ETFE of Y and U'}, 'Location', 'northwest');
xlabel('Frequency');
ylabel('Magnitude');
xlim([.2, max(omega(idx))]);
ylim([1e-2, 1e1]);
grid on

% Indirect approach (G*S and S) is worse than direct approach (Y and U)