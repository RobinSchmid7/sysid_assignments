%System Identification Ex 4
%Author: Robin Schmid, schmirob@ethz.ch
%% a) Generate y
close all; clear; clc; 
N = 1024;
t = 0:N-1;
e = sqrt(0.01)*randn(N,1); % Gaussian noise with variance 0.01

% Try different inputs
% omega_sin = 1;
% u = sin(omega_sin*2*pi/N.*x); % Option 1 for input: sinus
u = idinput(N); % Option 2 for input: random binary input
% u = randn(N,1); % Option 3 for input: normalized Gaussian
% u = u / max(abs(u));
%u = repmat(idinput(N/4),4,1) % Option 4 for input: periodic

% Simulate system
% Option 1 for tf
z = tf('z', 1); % T = 1 is the discrete sampling time, -1: for unspecified
G = 0.1*z/(z^4-2.2*z^3+2.42*z^2-1.87*z+0.7225);
H = 0.5*(z-0.9)/(z-0.25);
    % Option 2 for tf
    %G = tf([0,1, 0], [1, -2.2, 2.42 -1.87 0.7227], 1);

y_1 = lsim(G,u);
y_2 = lsim(H,e);
y = y_1 + y_2;

% Plot noise and response
figure(1)
plot(t, y);
hold on
plot(t, u);
xlim([0,N-1]);
legend({'y', 'u'});
grid on;

%% b) Estimate G and compare
% Plot G and H
Y = fft(y);
U = fft(u);

    % These spacings are equivalent
    % to not include the endpoint
    % o_1 = linspace(0,2*pi,N+1);
    % o_1 = o_1(1:end-1)';
    % o_2 = (2*pi/N)*[0:N-1]'

omega = (2*pi/N)*[0:N-1]';
idx = find(omega > 0 & omega < pi); % Only check for pos freq
    % Other option:
    % Y_N = fft(y); U_N = fft(u); G_N = Y_N(1:N/2+1)./U_N(1:N/2+1);
    % omega = 2*pi/N*(0:N/2); % Only positive frequencies
% loglog(omega(idx), abs(U(idx)));
% loglog(omega(idx), abs(Y(idx)));
G_est = Y./U;

% Freq response for certain omega
    % Option 1:
    % z = exp(1j.*omega);
    % P_freq = 0.1*z/(z^4-2.2*z^3+2.42*z^2-1.87*z+0.7225);
% Option 2:
G_freqresp = squeeze(freqresp(G, omega)); % Squeeze only extract the value for each freq
H_freqresp = squeeze(freqresp(H, omega));

% Plotting estimated response
figure(2)
tiledlayout(2,1);
nexttile;
loglog(omega(idx), abs(G_est(idx)), 'r-'); % Magnitude estimated
hold on % For plotting loglog use hold on after first plot!
loglog(omega(idx), abs(G_freqresp(idx)), 'b-'); % Magnitude true plant
loglog(omega(idx), abs(H_freqresp(idx)), 'g--'); % Magnitude true plant
xlim([.2,max(omega(N/2+1))]); % X limit not from 0 because of loglog plot
ylim([1e-2,1e2]);
legend({'ETFE unsmoothed', 'true plant', 'noise'});
xlabel('Frequency');
ylabel('Magnitude');
    % semilogx(omega(idx), angle(G_est(idx))); % For plotting phase

nexttile;
Error = G_est-G_freqresp;
loglog(omega(idx), abs(Error(idx))); % Magnitude error
xlim([.2,max(omega(N/2+1))]);
xlabel('Frequency');
ylabel('Magnitude');
title('Error response');

%% c) Split up data
N_batch = N/4;
u_reshaped = reshape(u, N_batch, 4); % Specify dim of new matrix, use [] to
% automatically calculate dim of new array
% u_respaped = [[u(1) u(5) u(9) u(13)]; [u(2) u(6) ...]; ...]; % Filled up
% row-wise
y_reshaped = reshape(y, N_batch, 4);

U_N_reshaped = fft(u_reshaped, [], 1); % Specify axis of fft, here along columns
Y_N_reshaped = fft(y_reshaped, [], 1);

G_est_split = Y_N_reshaped(1:N_batch/2+1, :)./U_N_reshaped(1:N_batch/2+1, :); % ETFE on pos freq for each batch
G_est = mean(G_est_split, 2); % Find average along all columns, i.e. along all batches
omega_split = 2*pi/N_batch*[0:N_batch/2]; % Frequencies of batch

% Plotting estimated response
figure(3)
tiledlayout(2,1);
nexttile;
loglog(omega_split, abs(G_est_split), 'r-'); % Magnitude estimated
hold on % For plotting loglog use hold on after first plot!
loglog(omega(idx), abs(G_freqresp(idx)), 'b-'); % Magnitude true plant
loglog(omega(idx), abs(H_freqresp(idx)), 'g--'); % Magnitude true plant
xlim([.2,max(omega(N/2+1))]); % X limit not from 0 because of loglog plot
ylim([1e-2,1e2]);
legend({'ETFE smoothed', 'true plant', 'noise'});
xlabel('Frequency');
ylabel('Magnitude');
    % semilogx(omega(idx), angle(G_est(idx))); % For plotting phase

nexttile;
G_freqresp = squeeze(freqresp(G, omega_split));
Error = G_est_split-G_freqresp;
loglog(omega(idx), abs(Error(idx))); % Magnitude error
xlim([.2,max(omega(N/2+1))]);
xlabel('Frequency');
ylabel('Magnitude');
title('Error response');