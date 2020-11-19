%-------------------------------------------------
% System Identification Exercise Set 7, Task 3
% Flavio Regenass
% November 2020
%-------------------------------------------------
% ARX Model

clc; clear all; close all;

% Parameters
a = 0.5;
b = 1;

var_1 = 0.2;
var_2 = 0.2;

%% Task a)

% Signal length
N_vec = [50, 250, 3000];
% N_vec = [10, 20, 50, 100, 200, 500, 1000, 2000];

% Number of Experiments
M = 1500;

% Tensor to store estimated parameters for M Experiments for each N-Value.
estimate_vector_all_N_1 = zeros(M,2,length(N_vec));
estimate_vector_all_N_2 = zeros(M,2,length(N_vec));

% Sum of squared Residuals, error of estimation
sum_residuals_1 = zeros(M, 1, length(N_vec));
sum_residuals_2 = zeros(M, 1, length(N_vec));

% Outer loop to go through Different Signal lengths
for k = N_vec
    % Input u
    u = randn(k,1);
    
    % Inner loop to conduct M experiments
    for j = 1:M
        
        % Noise realisations
        % a) Noise Normally distributed N(0, 0.2)
        w_1 = sqrt(var_1).*randn(k,1);
        % b) Noise uniformly distributed Uni(0, 0.2)
        w_2 = sqrt(var_2).*rand(k,1);
        
        % System Output with different noises
        y_1 = zeros(k,1);
        y_2 = zeros(k,1);
        for i = 2:k
            y_1(i) = a * y_1(i-1) + b * u(i-1) + w_1(i);
            y_2(i) = a * y_2(i-1) + b * u(i-1) + w_2(i);
        end
        
        % Regressor Matrix Phi
        phi_1 = zeros(k, 2);
        phi_2 = zeros(k, 2);
        for i = 2:k % Start at 2, since we assume the system to be at rest before we start (y = 0, u = 0)
            phi_1(i,:) = [y_1(i-1), u(i-1)];
            phi_2(i,:) = [y_2(i-1), u(i-1)];
        end
        
        % Least squares
        theta_hat_1 = phi_1 \ y_1;
        theta_hat_2 = phi_2 \ y_2;
        
        % Compute Residuals
        residuals_1 = (y_1 - phi_1 * theta_hat_1);
        residuals_2 = (y_2 - phi_2 * theta_hat_2);
        % Square them
        sq_residuals_1 = residuals_1.^2;
        sq_residuals_2 = residuals_2.^2;
        % Sum them up
        sum_residuals_1(j, 1, find(N_vec == k)) = 1/var_1 .* sum(sq_residuals_1);
        sum_residuals_2(j, 1, find(N_vec == k)) = 1/var_2 .* sum(sq_residuals_2);
        
        % Extend Estimate Vector
        estimate_vector_all_N_1(j, :, find(N_vec == k)) = [theta_hat_1(1), theta_hat_1(2)];

        estimate_vector_all_N_2(j, :, find(N_vec == k)) = [theta_hat_2(1), theta_hat_2(2)];
        
    end
end

%% Plots 
% Histogramme
% Case a) Normally Distributed Noise
% ========================================
figure;
for i = 1:length(N_vec)
    h_1_a = histogram(estimate_vector_all_N_1(:, 1, i));
%     legend(['N = ', num2str(N_vec(i))]);
    hold on
end
legend(['N = ', num2str(N_vec(1))], ['N = ', num2str(N_vec(2))], ['N = ', num2str(N_vec(3))]);
title('Estimates from M experiments for "a" with Data Length N, normally distributed noise');

figure;
for i = 1:length(N_vec)
    h_1_b = histogram(estimate_vector_all_N_1(:, 2, i));
%     legend(['N = ', num2str(N_vec(i))]);
    hold on
end
legend(['N = ', num2str(N_vec(1))], ['N = ', num2str(N_vec(2))], ['N = ', num2str(N_vec(3))]);
title('Estimates from M experiments for "b" with Data Length N, normally distributed noise')

% Case b) Uniformally distributed Noise
% =======================================
% Due to the nature of the uniformally distributed noise, there is no
% improvement in the estimation for longer N, since the influence of the
% noise does not decay..
% =======================================
figure;
for i = 1:length(N_vec)
    h_2_a = histogram(estimate_vector_all_N_2(:, 1, i));
%     legend(['N = ', num2str(N_vec(i))]);
    hold on
end
legend(['N = ', num2str(N_vec(1))], ['N = ', num2str(N_vec(2))], ['N = ', num2str(N_vec(3))]);
title('Estimates from M experiments for "a" with Data Length N, uniformally distributed noise');

figure;
for i = 1:length(N_vec)
    h_2_b = histogram(estimate_vector_all_N_2(:, 2, i));
%     legend(['N = ', num2str(N_vec(i))]);
    hold on
end
legend(['N = ', num2str(N_vec(1))], ['N = ', num2str(N_vec(2))], ['N = ', num2str(N_vec(3))]);
title('Estimates from M experiments for "b" with Data Length N, uniformally distributed noise')


%% Task b) Distribution of Residuals for Normally Distributed Noise

% Compute Squared Residuals
% --> See for loop above


% Normal Distribution Approximation
N = N_vec(end);
approx_normal = sqrt(2*(N-length(theta_hat_1))) .* randn(N, 1) + (N-length(theta_hat_1));

approx_chi2 = chi2pdf(linspace(0, N+1000),N-length(theta_hat_1));

figure;
h_sq_residuals_highN = histogram(sum_residuals_1(:,1, end)); % Histogram of Residuals
hold on
yyaxis left
h_approx = histogram(approx_normal); % approximate Chi Squared with Normal Distr.
hold on
yyaxis right
plot(linspace(0, N+1000), approx_chi2); % Chi Squared

legend('Distribution of Residuals for Estimate', 'Approximate Distribution from randn((N-p), 2(N-p))', 'PDF of Chi^2(N-p) distribution')
title(['Distribution of Residuals for N = ', num2str(N), ' Data Points and M = ', num2str(M), ' Experiments using Normally Distributed Noise']);
