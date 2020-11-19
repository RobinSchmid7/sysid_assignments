%-------------------------------------------------
% System Identification Exercise Set 2, Task 4
% Flavio Regenass
% October 2020
%-------------------------------------------------
% Estimate Multi-Spectral Transmissivity for a camera dependant on the
% wavelength of the absorbed radiation, given Camera-Measurements and other
% relevant data

%% Task 4, Part A
%----------------------------------------------------------------------------------------------------------------

%% Load the data
load('experiment3.mat');

%% Variables used
% Estimated Parameter Set theta_hat = [t(lambda_1), t(lambda_2), ..., t(lambda_N)]'

% Measurement Vector Y = [Y_R(T_i), Y_G(T_i), Y_B(T_i)]'
Y_R = table2array(Camera_measurements(1:end, 2));
Y_G = table2array(Camera_measurements(1:end, 3));
Y_B = table2array(Camera_measurements(1:end, 4));
Y = [Y_R; Y_G; Y_B];

% Temperature T
T = table2array(Camera_measurements(1:end, 1));
% Number of Temperatures:
n = length(table2array(Camera_measurements(1:end, 1)));

%Wavelength lambda
lambda = table2array(Camera_sensitivity(1:end, 1));
% Number of Wavelengths:
k = length(table2array(Camera_sensitivity(1:end, 1)));


% Camera Sensitivity s(lambda)
s_R = table2array(Camera_sensitivity(1:end, 2));
s_G = table2array(Camera_sensitivity(1:end, 3));
s_B = table2array(Camera_sensitivity(1:end, 4));

% Blackbody Emissivity e(T, Lambda)
emiss = table2array(BlackBody_emission(1:end, 2:end));

% Regressor PHI 
%(since we stack Y_R, Y_G, Y_B, we must take 3x theta_hat as
% parameter vector, since mxn = length(Y) x 3*16)
phi_R = zeros(length(Y_R), 16);
phi_G = zeros(length(Y_G), 16);
phi_B = zeros(length(Y_B), 16);

for i = 1:n
    for j = 1:k
    phi_R(i, j) = s_R(j) * emiss(i, j);
    phi_G(i, j) = s_G(j) * emiss(i, j);
    phi_B(i, j) = s_B(j) * emiss(i, j);
    end
end

phi = [phi_R; phi_G; phi_B];

%% Compute Multi-Spectral Transmittance t(lambda) = theta_hat

theta_hat = phi \ Y;
theta_hat_analytical = inv(phi'*phi)*phi'*Y;

%% Analyse Analytical Solution of Linear Least Squares Problem theta = inv(phi'*phi)*phi'*Y
% Compute Estimation with bad Theta
Y_hatAnalytical_R = phi_R* theta_hat_analytical;
Y_hatAnalytical_G = phi_G* theta_hat_analytical;
Y_hatAnalytical_B = phi_B* theta_hat_analytical;

%Compare Bad Estimate to true Camera Measurement
figure;
hold on
scatter(T, Y_R)
scatter(T, Y_G)
scatter(T, Y_B)
scatter(T, Y_hatAnalytical_R)
scatter(T, Y_hatAnalytical_G)
scatter(T, Y_hatAnalytical_B)
title('Comparison of Measurements Y and Estimates Y_{hatAnalytical}')
legend('Y_R', 'Y_G', 'Y_B', 'Y_{hat}_R', 'Y_{hat}_G', 'Y_{hat}_B');
xlabel('Temperature T');
hold off

%% Analyse Numerical Solution with Backslash Operator theta = phi \ Y 
% Compute Estimation with good Theta
Y_hat_R = phi_R* theta_hat;
Y_hat_G = phi_G* theta_hat;
Y_hat_B = phi_B* theta_hat;

% Compare good Estimate to true Camera Measurement
figure;
hold on
scatter(T, Y_R)
scatter(T, Y_G)
scatter(T, Y_B)
scatter(T, Y_hat_R)
scatter(T, Y_hat_G)
scatter(T, Y_hat_B)
title('Comparison of Measurements Y and Estimates Y_{hat}')
legend('Y_R', 'Y_G', 'Y_B', 'Y_{hat}_R', 'Y_{hat}_G', 'Y_{hat}_B');
xlabel('Temperature T');
hold off

figure;
scatter(lambda, theta_hat)
xlabel('Wavelength Lambda');
ylabel('Multi-Spectral Transmittance t(lambda)');
title('Unweighted, unregularized Estimate for Theta with normal least squares --> not good')

% No good fit to data? Looks good to  me...

%% Try Regularized Least Squares ( Minimise 1/2 *  norm(phi * theta - Y)^2 + lambda/2 * norm(theta)^2 with lambda being the regularization parameter)
% The analytical solution then becomes theta = inv(phi' * phi + lambda *
% I) * phi' * Y
regParamLambda = 0.9;
theta_hat_regularized = inv(phi' * phi + regParamLambda * eye(16)) * phi' * Y;

% Analyse Regularized Solution
% Compute Estimation with regularized Theta
Y_hatRegularized_R = phi_R* theta_hat_regularized;
Y_hatRegularized_G = phi_G* theta_hat_regularized;
Y_hatRegularized_B = phi_B* theta_hat_regularized;

% Compare good Estimate to true Camera Measurement
figure;
hold on
scatter(T, Y_R)
scatter(T, Y_G)
scatter(T, Y_B)
scatter(T, Y_hatRegularized_R)
scatter(T, Y_hatRegularized_G)
scatter(T, Y_hatRegularized_B)
title('Comparison of Measurements Y and Estimates Y_{hatRegularized}')
legend('Y_R', 'Y_G', 'Y_B', 'Y_{hatRegularized}_R', 'Y_{hatRegularized}_G', 'Y_{hatRegularized}_B');
xlabel('Temperature T');
hold off

figure;
scatter(lambda, theta_hat_regularized)
xlabel('Wavelength Lambda');
ylabel('Multi-Spectral Transmittance t(lambda)');
title('Solution with Regularized Least Squares')

%% Analyse Numerical Solution of Restricted Linear Least Squares Problem using lsqlin function
% theta = lsqlin(C, d, A, b, Aeq, beq, lb, ub);

C = phi;
d = Y;

lb = zeros(16, 1);
ub = ones(16, 1);

theta_hat_restricted = lsqlin(C, d, [], [], [], [], lb, ub);

% Compute Estimation with restricted Theta
Y_hatRestricted_R = phi_R* theta_hat_restricted;
Y_hatRestricted_G = phi_G* theta_hat_restricted;
Y_hatRestricted_B = phi_B* theta_hat_restricted;

%Compare Restricted Estimate to true Camera Measurement
figure;
hold on
scatter(T, Y_R)
scatter(T, Y_G)
scatter(T, Y_B)
scatter(T, Y_hatRestricted_R)
scatter(T, Y_hatRestricted_G)
scatter(T, Y_hatRestricted_B)
title('Comparison of Measurements Y and Estimates Y_{hatRestricted}')
legend('Y_R', 'Y_G', 'Y_B', 'Y_{hat}_R', 'Y_{hat}_G', 'Y_{hat}_B');
xlabel('Temperature T');
hold off

figure;
scatter(lambda, theta_hat_restricted);
xlabel('Wavelength Lambda');
ylabel('Multi-Spectral Transmittance t(lambda)');
title('Solution with Restricted Least Squares')

%% Compare different thetas
figure;
hold on
scatter(lambda, theta_hat);
scatter(lambda, theta_hat_analytical);
scatter(lambda, theta_hat_regularized);
scatter(lambda, theta_hat_restricted);
title('Comparison of Theta for different calculation Methods');
legend('theta_{hatNumerical}', 'theta_{hatAnalytical}', 'theta_{hatRegularized}', 'theta_{hatRestricted}');
hold off



