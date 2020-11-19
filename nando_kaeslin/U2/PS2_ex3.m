
%% Problem 3.1

A = importdata('experiment1.dat');

t = A.data(:,1);  % time
y = A.data(:,2);  % position
stddev = A.data(:,3);

close all
errorbar(t,y,stddev, 'bo')

k = length(t);

X = [ones(k,1) t 0.5*t.^2];

theta_hat = X\y;

t_LS = min(t):0.01:max(t);
y_LS = theta_hat(1) + theta_hat(2)*t_LS + 0.5*theta_hat(3)*t_LS.^2;

hold on
grid on
plot(t_LS, y_LS)
xlabel('t')
ylabel('y')


%% Problem 3.2

R = diag(stddev.^2);

cov_th_hat = (X'*X)\X'*R*X/(X'*X);

cov_th_hat1 = inv(X'/R*X); % this is the same as above

stddev_estimate = cov_th_hat.^.5;


%% Problem 3.3, 3.4

B = importdata('experiment2.dat');

stddev = B.data(:,3);

R = diag(stddev.^2);

W = inv(R);

theta_hat_weighted = (X'*W*X)\X'*W*y;

cov_th_hat_weighted = inv(X'/R*X);

stddev_estimate = cov_th_hat_weighted.^.5;
