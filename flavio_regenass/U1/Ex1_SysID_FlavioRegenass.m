% System Identification Exercise Set 1, Task 3
% Estimating theta from y = theta * x + 2 * theta * x^2 - theta * x^3 +
% v(k)
% Where v(k) is normally distributed with mu = 1.5 and sigma^2 = 1.2
% x and y measurements are given in "SysID_Exercise_1.mat"

%% Task 3, Part A
%------------------------------------------------------

%% Load the data
load('SysID_Exercise_1.mat');

%% Create Variables used
%System Input x = u
Z = y;
X = u;

mu = 1.5;
sigma = 1.2;

syms theta x z

%% Y Function to be analysed
yFunc = theta * (x + 2 * x^2 - x^3);

%% Compute Likelihood function
%likelihood = mle(Z, 'distribution', 'norm');

% Likelihood Function for sample z given theta, z being an element of total Range of samples Z
L = 1/(sqrt(2*pi*sigma)) * exp(-1/(2*sigma)*(z -(yFunc + mu))^2);

% Likelihood function to observe the whole of the samples Z given theta
totL = 1;
for i = 1:length(X)
    totL =totL * subs(L, [x, z], [X(i), Z(i)]);
end

% Logarithmical Likelihood function of totalLikelihood
logTotL = log(totL);

%% Derive log Likelihood func with respect to theta and set to zero to find best theta
diffLogL = diff(logTotL, theta);
% for n = 1:length(u)
%    likelihood = pdf(mu, sigma, Z);
% end

% Look for Theta maximising the diffLogL function
opt_theta = vpasolve(diffLogL == 0, theta)

%% Calculate resulting Mean Squared error when using theta
MAPerrSum = 0;
for i = 1:length(Z)
   MAPerrSum = MAPerrSum +(Z(i) - subs(yFunc, [x, theta], [X(i), opt_theta]))^2;
end
MSE_withoutMAP = 1/length(Z) * MAPerrSum;


%% Task 3, Part B
% Maximum a Posteriori Estimate for Theta
%------------------------------------------------------

%% A priori Parameter Distribution

% Parameters for Parameter Distribution
mu_theta = 1.5;
sigma_theta = 0.3^2;

% New Likelihood function MAPtotL
% L = totL * f_theta(theta) with theta being the sought parameter
% With f_theta = (1/(sqrt(2*pi*sigma_theta)) *
% exp(-(1/(2*sigma_theta))*(theta - mu_sigma)^2) --> Slide 2.22 Lecture 2

f_theta = 1/(sqrt(2*pi*sigma_theta)) * exp(-(1/(2*sigma_theta))*(theta - mu_theta)^2);
MAPtotL = totL * f_theta;

% Build Log Likelihood function
% Logarithmical Likelihood function of totalLikelihood
MAPlogTotL = log(MAPtotL);

%% Derive log Likelihood func with respect to theta and set to zero to find best theta
diffMAPLogL = diff(MAPlogTotL, theta);
% for n = 1:length(u)
%    likelihood = pdf(mu, sigma, Z);
% end

% Look for Theta maximising the diffLogL function
MAP_opt_theta = vpasolve(diffMAPLogL == 0, theta)

%% Calculate resulting Mean Squared error when using MAP theta
MAPerrSum = 0;
for i = 1:length(Z)
   MAPerrSum = MAPerrSum +(Z(i) - subs(yFunc, [x, theta], [X(i), MAP_opt_theta]))^2;
end

%Compare the errors of the two estimations
MSE_withoutMAP
MSE_withMAP = 1/length(Z) * MAPerrSum




