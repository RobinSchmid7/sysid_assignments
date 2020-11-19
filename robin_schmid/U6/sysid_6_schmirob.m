%System Identification Ex 6
%Author: Robin Schmid, schmirob@ethz.ch
%% 1. Apply 2 periods and average
clear all; close all; clc;
tau_max = 80;
K = 2*tau_max; % Calculation length

% Prbs signal with length 80, magnitude 2
% Attention: prbs function can only create signals with length 2^n-1
% This if the length is not equal to 2^n-1 then the prbs signal is cut and
% the resulting autocorrelation is not periodic
u_period1 = 2*idinput(tau_max,'PRBS');
% Thus better use a random binary signal in the range and scale it
% accordingly, guaranteed to have a period of tau_max
u_period2 = (4*randi([1,2],1,tau_max)-6)';
% plot(1:size(u_period2,2),u_period2);

% Compare effect on autocorrelation of both options
R_u1 = corr_per(u_period1, u_period1);
R_u2 = corr_per(u_period2, u_period2);

f1 = figure(1);
set(f1, 'visible', 'off');
plot(1:80, R_u1, 'b-');
hold on
plot(1:80, R_u2, 'r-');
title('Autocorrelation for 1: idinput() and 2: randi() approach');
xlabel('Time');
ylabel('u');
legend({'R_{u1}', 'R_{u2}'});
grid on
% Here both methods give similar results
% Use second approach for u_period

u_period = u_period2;
u = [u_period; u_period];
% plot(1:size(u,1),u);

% Create Toeplitz matrix
Phi_u = toeplitz(u, zeros(1,tau_max)); % Usage: columns, rows

% Two experiments
g_est = zeros(tau_max,2);
for i = 1:2
    % Noise distrubution, assume Gaussian with variance 0.05
    v = sqrt(0.05)*randn(K, 1);

    y = zeros(K,1);
    y(1) = v(1);
    for t = 1:K-1
        y(t+1) = 5/8*y(t) + 11/10*u(t) + v(t+1);
    end
    % plot(1:size(y,1),y)
    
    g_est(:,i) = Phi_u\y;
end

% Average estimates
g_mean = mean(g_est,2);

%% 2. Apply N periods
tau_max = 80;
u_period = 2*idinput(tau_max);

% % Noiseless pulse response, unstable
% y = zeros(tau_max,1);
% for t = 1:tau_max-1
%     y(t+1) = 5/8*y(t) + 11/10*u(t);
% end
% Phi_u_0 = toeplitz(u_period, zeros(1,tau_max))
% g_est_0 = Phi_u_0\y

% For N periods
g_est_N = zeros(tau_max,2);
err_N = zeros(1,30);

for i = 2:30
    K = i*tau_max;
    v = sqrt(0.05)*randn(K, 1);
    u = repmat(u_period, 1, i);
    
    Phi_u = toeplitz(u, zeros(1,tau_max));
    
    y_0 = zeros(K,1);
    y_0(1) = 0;
    for t = 1:K-1
        y_0(t+1) = 5/8*y_0(t) + 11/10*u(t);
    end
    
    y = zeros(K,1);
    y(1) = v(1);
    for t = 1:K-1
        y(t+1) = 5/8*y(t) + 11/10*u(t) + v(t+1);
    end
    
    g_est_0 = Phi_u\y_0;
    g_est_N(:,i) = Phi_u\y;
    
    % Plotting
    err_N1(i) = norm(g_est_0 - g_est_N(:,i))^2;
end

% Error for averaging of part 1
err_avg = norm(g_est_0(1:tau_max) - g_mean)^2;

f2 = figure(2);
set(f2, 'visible', 'off');
plot(1:30, err_N1, 'b-');
hold on
yline(err_avg, 'r-');
title('Error for N periods with prbs input');
xlabel('Periods');
ylabel('Error');
grid on

% Is sentive to noise, if converging at around 5-10 period method 2 is better
% than method 1

%% 3. Use uniform random noise as input
tau_max = 80;
K = 2*tau_max; % Calculation length

% Prbs signal with length 80, magnitude 2
u_period = 2*rand(tau_max,1);
u = [u_period; u_period];
% plot(1:size(u,1),u)

% Create Toeplitz matrix
Phi_u = toeplitz(u, zeros(1,tau_max)); % Usage: columns, rows

% Two experiments
g_est = zeros(tau_max,2);
for i = 1:2
    % Noise distrubution, assume Gaussian with variance 0.05
    v = sqrt(0.05)*randn(K, 1);

    y = zeros(K,1);
    for t = 1:K-1
        y(t+1) = 5/8*y(t) + 11/10*u(t) + v(t+1);
    end
    % plot(1:size(y,1),y)
    
    g_est(:,i) = Phi_u\y;
end

% Average estimates
g_avg = mean(g_est,2);

%
tau_max = 80;
u_period = 2*idinput(tau_max);

% % Noiseless pulse response, unstable
% y = zeros(tau_max,1);
% for t = 1:tau_max-1
%     y(t+1) = 5/8*y(t) + 11/10*u(t);
% end
% Phi_u_0 = toeplitz(u_period, zeros(1,tau_max))
% g_est_0 = Phi_u_0\y

% For N periods
g_est_N = zeros(tau_max,2);
err_N = zeros(1,30);

for i = 2:30
    K = i*tau_max;
    v = sqrt(0.05)*randn(K, 1);
    u = repmat(u_period, 1, i);
    
    Phi_u = toeplitz(u, zeros(1,tau_max));
    
    y_0 = zeros(K,1);
    for t = 1:K-1
        y_0(t+1) = 5/8*y_0(t) + 11/10*u(t);
    end
    
    y = zeros(K,1);
    for t = 1:K-1
        y(t+1) = 5/8*y(t) + 11/10*u(t) + v(t+1);
    end
    
    g_est_0 = Phi_u\y_0;
    g_est_N(:,i) = Phi_u\y;
    
    % Plotting
    err_N2(i) = norm(g_est_0 - g_est_N(:,i))^2;
end

% Error for averaging of part 1
err_avg = norm(g_est_0(1:tau_max) - g_mean)^2;

% Compare all three methods
f3 = figure(3);
set(f3, 'visible', 'on');
p1 = plot(1:30, err_N1, 'b-');
hold on
p2 = yline(err_avg, 'y-');
p3 = plot(1:30, err_N2, 'r-');
title('Error for N periods with averaging and prbs and uniform noise');
xlabel('Periods');
ylabel('Error');
legend([p2, p1, p3], {'averaging', 'prbs', 'uniform'});
grid on
hold off

% Random input worse than using a periodic prbs signal

% Autocorrelation function to see effect of truncation on u_period
function [R] = corr_per(x, u)
    N = size(u, 1);
    R = zeros(1, N);
    
    % Check if input signal is even
    if mod(N, 2) == 0
        lags = -N/2+1:N/2;
    else
        % Shift such that lags are integer numbers
        lags = -(N-1)/2:(N-1)/2;
    end
    
    for t = 1:N
        for k = 1:N
            % Account for periodicity of signal
            if k-lags(t)<=0
                R(t) = R(t) + 1/N * x(k) * u(k-lags(t)+N);
            elseif k-lags(t)>N
                R(t) = R(t) + 1/N * x(k) * u(k-lags(t)-N);
            else
                R(t) = R(t) + 1/N * x(k) * u(k-lags(t));
            end
        end
    end
end
