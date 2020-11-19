%System Identification Ex 2
%Author: Robin Schmid, schmirob@ethz.ch
%% 3 (maximum likelihood, uncorrelated noise)
% Intitialize
importdata('experiment1.dat');
t = experiment1.t;
y = experiment1.y;
stdy = experiment1.stdy;
A = [1*ones(10,1) t 0.5*t.^2];

% Solve unweighted least squares problem
p_est = A\y;
y_est = p_est(1) + p_est(2)*t_step + p_est(3)*0.5*t_step.^2;

% Plotting
figure(1)
errorbar(t, y, stdy, 'bo') % plot with variance
hold on
t_step = min(t):0.01:max(t); % use finer sampling for time
plot(t_step, y_est)
legend('measured', 'estimated')

% Find estimated mean and variance
cov = inv(A'*A)*stdy(1)^2;
est_mean = [p_est(1); p_est(2); p_est(3)]
est_var = [sqrt(cov(1,1)); sqrt(cov(2,2)); sqrt(cov(3,3))]
%% 3 (maximum likelihood, correlated noise)
% Initialize
importdata('experiment2.dat');
t = experiment2.t;
y = experiment2.y;
stdy = experiment2.stdy;
R = diag(stdy.^2);
P = [1*ones(10,1) t 0.5*t.^2];

% Solve weighted least square problem
A = P'*inv(R)*P;
b = P'*inv(R)*y;
p_est = A\b;

% Find estimated mean and variance
cov = inv(P'*inv(R)*P);
est_mean = [p_est(1); p_est(2); p_est(3)]
est_var = [sqrt(cov(1,1)); sqrt(cov(2,2)); sqrt(cov(3,3))]
%% 4 (regularized least squares)
% Initialize
load('experiment3.mat');
s_r = Camera_sensitivity.Red;
s_g = Camera_sensitivity.Green;
s_b = Camera_sensitivity.Blue;
lambda = Camera_sensitivity.Wavelenght;
e = BlackBody_emission{:,2:17};
y_r = Camera_measurements.Red;
y_g = Camera_measurements.Green;
y_b = Camera_measurements.Blue;

% Create matrices for each color
A_r = zeros(151,16);
for i = 1:size(e,1)
    for j = 1:size(e,2)
    A_r(i,j) = s_r(j)*e(i,j);
    end
end

A_g = zeros(151,16);
for i = 1:size(e,1)
    for j = 1:size(e,2)
    A_g(i,j) = s_g(j)*e(i,j);
    end
end

A_b = zeros(151,16);
for i = 1:size(e,1)
    for j = 1:size(e,2)
    A_b(i,j) = s_b(j)*e(i,j);
    end
end

% Do not use standart least squares, use constraint least square
% x_r = lsqlin(A_r,y_r,[],[],[],[],zeros(16,1),ones(16,1)); % Here solve individual least square problems
% x_b = lsqlin(A_b,y_b,[],[],[],[],zeros(16,1),ones(16,1));
% x_g = lsqlin(A_g,y_g,[],[],[],[],zeros(16,1),ones(16,1));

% Solve stacked least square problem
y_tot = [y_r; y_g; y_b];
A_tot = [A_r; A_g; A_b];
%x_tot = A_tot\y_tot; % If not using constraint problem get non-physical results
x_tot = lsqlin(A_tot,y_tot,[],[],[],[],zeros(16,1),ones(16,1));
plot(lambda, x_tot, '*-')
xlabel('Wavelength (nm)')
ylabel('Transmittance')