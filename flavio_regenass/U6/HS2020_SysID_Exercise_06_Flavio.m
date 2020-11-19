%-------------------------------------------------
% System Identification Exercise Set 6, Task 3
% Flavio Regenass
% November 2020
%-------------------------------------------------
% Pulse Response

clc; clear all;

M = 80; % Period Length of PRBS Signal
tau_max = 80;
m = 2; % number of Periods

% Parameters
a = 5/8;
b = 11/10;

%% Task 1: Averaging

% PRBS input
u_seq = 4*round(rand(1, M)) - 2; % 1 Period
u = [];
for k = 1:m
    u = [u, u_seq]; % Assemble Periodical Input
end

% Form Toeplitz Matrix
T = toeplitz(u);
T = T - triu(T, 1); % Subtract upper right part for causal System (all zeros)

% Truncate Toeplitz Matrix for Tau_max
phi = T(:, 1:tau_max+1);

g_hat_avg = [];
g_hat = zeros(tau_max+1, 2);

for i = 1:2 % Iterate twice and Average
    
    % Noise
    cov = 0.05 * eye(m*M);
    mean = zeros(m*M,1);
    v = mvnrnd(mean,cov,1);
    
    % Output y
    y = zeros(m*M,1);
    
    for t = 1:m*M
        if t == 1
            y(1) = a * 0 + b * 0 + v(1);
        else
            y(t) = a * y(t-1) + b * u(t-1) + v(t);
        end
    end
    % Estimate Transfer Function
    g_hat(:, i) = (phi.' * phi) \ phi.' * y;
end

g_hat_avg = 1/2 * (g_hat(:, 1) + g_hat(:, 2));

%% Task 2: Increase Measurement Length K

m_vec = [2:30]; % Number of Periods measured
MSE2 = zeros(1,length(m_vec));
g_hat_task2 = zeros(tau_max+1, length(m_vec));

% Compute true Transfer function g_0, without noise added
% Output y
y_0 = zeros(m*M,1);
for t = 1:m*M
    if t == 1
        y_0(1) = 0;
    else
        y_0(t) = a * y(t-1) + b * u(t-1);
    end
end
% Estimate Transfer Function
g_0 = (phi.' * phi) \ phi.' * y_0;

% Compute estimates for varying Period Counts Applied
for periodCount = m_vec
    
    % PRBS input
    u_seq = 4*round(rand(1, M)) - 2; % 1 Period
    u2 = [];
    for k = 1:periodCount
        u2 = [u2, u_seq]; % Assemble Periodical Input
    end
    
    % Form Toeplitz Matrix
    T2 = toeplitz(u2);
    T2 = T2 - triu(T2, 1); % Subtract upper right part for causal System (all zeros)
    
    % Truncate Toeplitz Matrix for Tau_max
    phi = T2(:, 1:tau_max+1);
    
    % Noise
    cov = 0.05 * eye(periodCount*M);
    mean = zeros(periodCount*M,1);
    v = mvnrnd(mean,cov,1);
    
    % Output y
    y2 = zeros(periodCount*M,1);
    
    for t = 1:periodCount*M
        if t == 1
            y2(1) = a * 0 + b * 0 + v(1);
        else
            y2(t) = a * y2(t-1) + b * u2(t-1) + v(t);
        end
    end
        
    % Estimate Transfer Function
    g_hat_task2(:, periodCount) = (phi.' * phi) \ phi.' * y2;
    
    % Calculate the error of the Estimate
    MSE2(periodCount-1) = sqrt(sum((g_hat_task2(:, periodCount) - g_0).^2));
end

% Plot the result
figure;
plot(m_vec, MSE2)
hold on
plot(2, sqrt(sum((g_hat_avg - g_0).^2)), 'r*')
legend('Estimates with different Period Counts', 'Averaged Estimate with K = 2 Periods')
xlabel('Number of Periods m in data K = m * M')
ylabel('MSE compared to TF without Noise (g_{0})')
title('Comparing MSE of different Period-Counts with Averaged Estimate')

%% Task 3: Make Signal longer

MSE3 = zeros(1,length(m_vec));
g_hat_task3 = zeros(tau_max+1, length(m_vec));

% Compute estimates for varying Period Counts Applied
for periodCount = m_vec
    
    % PRBS input
    u_seq3 = 4*round(rand(1, periodCount*tau_max)) - 2; % 1 Period
    u3 = u_seq3;
    
    % Form Toeplitz Matrix
    T3 = toeplitz(u3);
    T3 = T3 - triu(T3, 1); % Subtract upper right part for causal System (all zeros)
    
    % Truncate Toeplitz Matrix for Tau_max
    phi3 = T3(:, 1:tau_max+1);
    
    % Noise
    cov = 0.05 * eye(periodCount*M);
    mean = zeros(periodCount*M,1);
    v = mvnrnd(mean,cov,1);
    
    % Output y
    y3 = zeros(periodCount*tau_max,1);
    
    for t = 1:periodCount*M
        if t == 1
            y3(1) = a * 0 + b * 0 + v(1);
        else
            y3(t) = a * y3(t-1) + b * u3(t-1) + v(t);
        end
    end
        
    % Estimate Transfer Function
    g_hat_task3(:, periodCount) = (phi3.' * phi3) \ phi3.' * y3;
    
    % Calculate the error of the Estimate
    MSE3(periodCount-1) = sqrt(sum((g_hat_task3(:, periodCount) - g_0).^2));
end

% Plot the result
figure;
plot(m_vec, MSE3)
hold on
plot(2, sqrt(sum((g_hat_avg - g_0).^2)), 'r*')
legend('Estimates with different Signal Lengths', 'Averaged Estimate with K = 2 Periods')
xlabel('Data Length Index (N) K = N * tau_max')
ylabel('MSE compared to TF without Noise (g_{0})')
title('Comparing MSE of different Signal Lengths with Averaged Estimate')

