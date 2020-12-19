%System Identification Exercise Session 5
%Author: Robin Schmid, schmirob@ethz.ch
close all; clear; clc;
[u_prbs,y_prbs,u_rand,y_rand] = GenerateData();

N = size(u_prbs,1); % 1023

% PRBS
Y_prbs = fft(y_prbs);
U_prbs = fft(u_prbs);
G_prbs = Y_prbs./U_prbs; % Unsmoothed

% Rand
Y_rand = fft(y_rand);
U_rand = fft(u_rand);
G_rand = Y_rand./U_rand; % Unsmoothed

G_est = G_prbs;
U = U_prbs;
G_s = zeros(N,1);

% For prbs gamma = 100 is good
% For rand gamma = 60 is good
% Reason: gamma proportional to variance of noise, for prbs have a small
% SNR, thus the noise is relatively small and can use a big gamma until the
% effect of the noise becomes significant -> can use bigger smoothing than
% with a rand signal -> this is an advantage of a prbs signal
figure(1)
it = 0;
for gamma = 80:20:150
    it = it+1;
    % Smoothing
    [omega, Wg] = WfHann(gamma, N);
    zidx = find(omega==0);
    omega = [omega(zidx:N,1); omega(1:zidx-1,1)];
    Wg = [Wg(1,zidx:N) Wg(1,1:zidx-1)];
    a = U.*conj(U); % Magnitude of U squared, variance weigthing

    for wn = 1:N
        Wnorm = 0;
        for xi = 1:N
            widx = mod(xi-wn,N)+1;
            G_s(wn) = G_s(wn) + Wg(widx) * G_est(xi) * a(xi);
            Wnorm = Wnorm + Wg(widx) * a(xi);
        end
        G_s(wn) = G_s(wn)/Wnorm;
    end
    loglog(omega, abs(G_s)); % Magnitude estimated
    legendInfo{it} = ['gamma = ' num2str(gamma)];
    hold on
end

legendInfo{it+1} = 'umsmoothed';
% Plot unsmoothed
loglog(omega, abs(G_prbs));
legend(legendInfo(:));