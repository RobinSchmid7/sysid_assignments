% System Identification Exercise Session 7
% Robin Schmid, schmirob@eth.ch
%%
clear all; close all; clc;

a = 0.9;
b = 0.5;
N_exp = 5:100;

% ARX model
z = tf('z', 1);
A = 1-a*z^-1;
B = b*z^-1;

% Noise and input Gaussian
u_init = 0.25*randn(1,100);
v_init = 0.25*randn(1,100);

theta0 = zeros(2,100);
theta_n = zeros(2,100);

it = 0;
% Run for same noise
for N = N_exp
    u = u_init(1:N);
    v = v_init(1:N);

    % Without noise
    y0 = lsim(B/A,u);
    % With noise
    y_n = lsim(B/A,u) + lsim(1/A,v);

    % LS for noise free case, assume at rest
    Phi = zeros(N,2);
    Phi(1,:) = [0 0];
    for i = 2:N
        Phi(i,:) = [y0(i-1) u(i-1)];
    end
    
    theta0(:,N) = Phi\y0;
    
    % LS with noise, assume at rest
    Phi = zeros(N,2);
    Phi(1,:) = [0 0];
    for i = 2:N
        Phi(i,:) = [y_n(i-1) u(i-1)];
    end
    
    theta_n(:,N) = Phi\y_n;
end
theta0(:,1:4) = [];
theta_n(:,1:4) = [];

% Plotting
figure(1);
subplot(2,1,1);
plot(N_exp,theta_n(1,:),'r');
hold on
yline(a,'b');
title('Estimate of a for different N');
xlabel('N');
ylabel('a');

subplot(2,1,2);
plot(N_exp,theta_n(2,:),'r');
hold on
yline(b,'b');
title('Estimate of b for different N');
xlabel('N');
ylabel('b');

% %% For multiple N
% l_s = 5:100;
% num_N = size(l_s,2);
% 
% for k = 1:num_N
%     N = l_s(k);
%     u = randn(N,1)*sqrt(var);
%     % w = randn(N,1)*sqrt(0.2);
%     w = randn(N,1)*sqrt(var);
%     % w = zeros(N,1);
%     phi = zeros(N,2);
%     y = zeros(N,1);
%     
%     y(1) = w(1); % Initialize y with first value
%     for t = 2:N
%         y(t) = a*y(t-1) + b*u(t-1) + w(t);
%         phi(t,:) = [y(t-1), u(t-1)];
%     end
%     theta_hat(k, :) = (phi'*phi)\phi'*y;
% end
% 
% figure(1);
% subplot(2,1,1);
% plot(l_s,theta_hat(:,1),'r');
% hold on
% yline(a,'b');
% title('Estimate of a for different N');
% xlabel('N');
% ylabel('a');
% 
% subplot(2,1,2);
% plot(l_s,theta_hat(:,2),'r');
% hold on
% yline(b,'b');
% title('Estimate of b for different N');
% xlabel('N');
% ylabel('b');