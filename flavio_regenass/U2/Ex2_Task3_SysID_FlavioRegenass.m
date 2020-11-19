%-------------------------------------------------
% System Identification Exercise Set 2, Task 3
% Flavio Regenass
% October 2020
%-------------------------------------------------
% Gravitational Experiment in Space to define g-constant of Planet
% Linear least squares to approx. y(t) = y_0 + v_0 * t + 1/2 * a * t^2

%% Task 3, Part A
%----------------------------------------------------------------------------------------------------------------

%% Load the data
% I used the "Import Data" function manually to import "experiment1.dat"
% and "experiment2.dat"

%% Variables used
% Estimated Parameter Set theta_hat = [y_0, v_0, 1/2*a]'
% Measurement Vector Y = [y_0 ... y_N]'
Y = table2array(experiment1(1:end, 2));
% timesteps t_i
T = table2array(experiment1(1:end, 1));
% Regressor PHI
phi = zeros(10, 3);
phi(1:end, 1) = 1;
phi(1:end, 2) = T;
phi(1:end, 3) = T.^2;

%% Estimate Parameters
% As Recommended in Lecture: Theta = phi' \ Y
theta_hat_numerical = phi\Y;

% As constructed analytically theta_hat = (phi' * phi)^-1 * phi' * Y
theta_hat_analytical = (phi.'*phi) \ phi.'*Y;

% --> Both theta_hat are equal! It might make a difference for the computer
% in terms of computational effort.

% According to the Estimate we are on the moon G_R!

%% Estimate Y_hat(t) with calculated parameters
Y_hat = phi * theta_hat_numerical;

figure;
scatter(T, Y)
hold on
scatter(T, Y_hat)
legend('Measurements', 'Estimations with LSE-Params');
title('Experiment 1')
hold off

%% Task 3, Part B
%-------------------------------------------------------------------------------------------------------------------
% Calculate Standard Deviation of the estimates for theta_hat
% Given in Slides: cov(theta_hat) = phi' * E(residual) * phi, where the
% E(residual) is equal to the variance of the gravity-constant in this
% example
% Slide 3.23
R_1 = diag(table2array(experiment1(1:end, 3)).^2);
%pre_Term1 = ((phi.' * phi) \ phi.').';

%covariance = pre_Term1.' * R_1 * pre_Term1;
covariance = ((phi.' * phi) \ phi.') * R_1 * phi / (phi.' * phi);
covariance_BLUE = inv(phi.' * (R_1 \ phi))

% Standard Deviation
stdev_exp1 = [sqrt(covariance(1,1)), sqrt(covariance(2,2)), sqrt(covariance(3,3))]
stdev_BLUE_exp1 = sqrt(covariance_BLUE)

%% Task 3, Part C
%-------------------------------------------------------------------------------------------------------------------
% Same as above, but with new data and weighted least squares.

%% Variables used
% Estimated Parameter Set theta_hat = [y_0, v_0, 1/2*a]'
% Measurement Vector Y = [y_0 ... y_N]'
Y_exp2 = table2array(experiment2(1:end, 2));
% timesteps t_i
T_exp2 = table2array(experiment2(1:end, 1));
% Regressor PHI
phi_exp2 = zeros(10, 3);
phi_exp2(1:end, 1) = 1;
phi_exp2(1:end, 2) = T_exp2;
phi_exp2(1:end, 3) = T_exp2.^2;

% Extract std dev and compute weights
stdev_data_exp2 = table2array(experiment2(1:end, 3));
variance_exp2 = stdev_data_exp2.^2;

% Weight Matrix W is diagonal and consists of inverse of variance
W = diag(variance_exp2.^(-1));


%% Estimate Parameters

% As constructed analytically theta_hat = (phi' * phi)^-1 * phi' * Y
theta_hat_exp2_analytical = (phi_exp2.' * W * phi_exp2) \ phi_exp2.' * W * Y_exp2;

% According to the Estimate we are on the moon G_I!

%% Estimate Y_hat(t) with calculated parameters
Y_hat_exp2 = phi_exp2 * theta_hat_exp2_analytical;

figure;
scatter(T_exp2, Y_exp2)
hold on
scatter(T_exp2, Y_hat_exp2)
legend('Measurements', 'Estimations with LSE-Params');
title('Experiment 2')

%% Task 3, Part D
%-----------------------------------------------------------------------------------------------------------------
% Calculate Standard Deviation of the estimates for theta_hat_exp2
% Given in Slides: cov(theta_hat) = phi' * E(residual) * phi, where the
% E(residual) is equal to the variance of the gravity-constant in this
% example
% Slide 3.23

R_2 = diag(table2array(experiment2(1:end, 3)).^2);
pre_Term2 = (phi_exp2.' * phi_exp2) \ phi_exp2.';

covariance_exp2 = pre_Term2 * R_2 * pre_Term2.';

% Standard Deviation
stdev_exp2 = [sqrt(abs(covariance_exp2(1,1))), sqrt(abs(covariance_exp2(2,2))), sqrt(abs(covariance_exp2(3,3)))]

% --> Standard Deviation became much smaller. Why?



