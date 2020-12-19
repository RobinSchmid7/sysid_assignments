% System Identification Ex 8
% Robin Schmid, schmirob@ethz.ch
%% 1. + 2. + 3.
clear all; close all; clc;

% Initialization
N = 1e4;
t = 0:N-1;

E = 100; % Number of experiments

Ts = -1; % Use unspecified sampling time (use value?)
z = tf('z', Ts);

% Tf for process
A = 1-1.5*z^-1+0.7*z^-2;
B = z^-1+0.5*z^-2;
C = 1-z^-1+0.2*z^-2;

% Tf for input
A_u = 1+0.1*z^-1-0.12*z^-2;
B_u = z^-1+0.2*z^-2;

% For multiple e, also change e_u each time (assume yes)?
theta_est = zeros(E,4);
for k = 1:E
    
    e = randn(N,1);
    e_u = randn(N,1);

    % Generate input sequence
    u = lsim(B_u/A_u, e_u, t);
    % plot(1:size(u,1),u);
    % Generating output
    y_e_true = lsim(C/A, e, t);
    y_u_true = lsim(B/A, u, t);
    y_true = y_e_true + y_u_true;

    % plot(1:size(y_true,1),y_true,'-');

    % Assume filter C given, find y_f through filtering
    % Work with this data, for filtered system work with ARX model
    y_f = lsim(1/C, y_true, t);
    u_f = lsim(1/C, u, t);

    % Estimating A_est and B_est
    Phi = zeros(N,4);
    Phi(1,:) = [0 0 0 0]; % Assume system is at rest at beginning
    Phi(2,:) = [y_f(1) 0 u_f(1) 0];

    for i = 3:N
        Phi(i,:) = [-y_f(i-1) -y_f(i-2) u_f(i-1) u_f(i-2)];
    end

    theta_est(k,:) = Phi\y_f;

    % Estimates for filtered output
    A_est = 1+theta_est(k,1)*z^-1+theta_est(k,2)*z^-2;
    B_est = theta_est(k,3)*z^-1+theta_est(k,4)*z^-2;

    % Option 1: Generate estimated filtered output, system without noise
    y_est = lsim(B_est/A_est, u, t); 
    % Option 2: Use direct mapping of filtered output
    % y_f = Phi*theta_est;
    % y_est = lsim(C, y_f, t); 
end

% Plotting
figure(1);
plot(1:N, y_true, 'b-');
hold on
plot(1:N, y_est, 'r-');
legend({'True output', 'Estimated output'});
title('True and estimated output');
xlabel('Time');
ylabel('Magnitude');
% Very good fit with this method, use smaller N to better see the quality

% Histogram for multiple experiments
figure(2);
subplot(2,2,1);
histogram(theta_est(:,1));
title('a1');
subplot(2,2,2);
histogram(theta_est(:,2));
title('a2');
subplot(2,2,3);
histogram(theta_est(:,3));
title('b1');
subplot(2,2,4);
histogram(theta_est(:,4));
title('b2');
% Values seem to be unbiased