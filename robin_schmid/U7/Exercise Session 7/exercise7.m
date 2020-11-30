% System Identification Exercise Session 7
% Robin Schmid, schmirob@eth.ch
%%
clear all; close all; clc;
a = 0.9;
b = 0.5;
var = 1;
N = 20;
%% For multiple N
l_s = 5:100;
num_N = size(l_s,2);

for k = 1:num_N
    N = l_s(k);
    u = randn(N,1)*sqrt(var);
    % w = randn(N,1)*sqrt(0.2);
    w = randn(N,1)*sqrt(var);
    % w = zeros(N,1);
    phi = zeros(N,2);
    y = zeros(N,1);
    
    y(1) = w(1); % Initialize y with first value
    for t = 2:N
        y(t) = a*y(t-1) + b*u(t-1) + w(t);
        phi(t,:) = [y(t-1), u(t-1)];
    end
    theta_hat(k, :) = (phi'*phi)\phi'*y;
end

figure(1);
subplot(2,1,1);
plot(l_s,theta_hat(:,1),'r');
hold on
yline(a,'b');
title('Estimate of a for different N');
xlabel('N');
ylabel('a');

subplot(2,1,2);
plot(l_s,theta_hat(:,2),'r');
hold on
yline(b,'b');
title('Estimate of b for different N');
xlabel('N');
ylabel('b');




