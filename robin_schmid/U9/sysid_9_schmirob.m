% System Identification Ex 8
% Robin Schmid, schmirob@ethz.ch
%%
clear all; close all; clc;

% 1. Generate data
N = 200;
t = 0:N-1;

% Use a different seed to generate input and error noise s.t. they are
% different
u = randn(1,N);

threshold = 1000;

Ts = 0.1;
z = tf('z', Ts);

A_true = 1 - 1.1*z^-1 + 0.18*z^-2;
B_true = z^-1 - 2*z^-2;

y_true = lsim(B_true/A_true, u); % No noise for true system?

models = [1, 2; 2, 2; 3, 2; 2, 1; 4, 3];
n_models = size(models,1);

MSE_LS = zeros(1,n_models);
MSE_IV = zeros(1,n_models);

y_est_LS = zeros(N,n_models);
y_est_IV = zeros(N,n_models);

for k = 1:n_models
    n = models(k,1);
    m = models(k,2);
    
    % 2. LS method (biased)
    e_LS = sqrt(0.5)*randn(1,N);
    y_LS = lsim(B_true/A_true, u) + e_LS'; % Attention lsim returns column vector
    
    % Phiy
    Phiy = zeros(N,n);    
    for iy = 1:N
        for jy = 1:n
            % Use at rest assumption
            if (iy-jy) < 1
                Phiy(iy,jy) = 0;
            else
                Phiy(iy,jy) = -y_LS(iy-jy);
            end
        end
    end
    
    % Phiu
    Phiu = zeros(N,m);
    for iu = 1:N
        for ju = 1:m
            % Use at rest assumption
            if (iu-ju) < 1
                Phiu(iu,ju) = 0;
            else
                Phiu(iu,ju) = u(iu-ju);
            end
        end
    end
    
    Phi = [Phiy Phiu];
    % Least square regression
    theta_LS = Phi\y_LS;
    
    % A_LS
    A_LS = 1;
    for iy = 1:n
        A_LS = A_LS + theta_LS(iy) * z^(-iy);
    end
    
    % B_LS
    B_LS = 0;
    for iu = 1:m
        B_LS = B_LS + theta_LS(n+iu) * z^(-iu);
    end
    
    % 3. Generate new data, use same u
    % Generate instrumental variable
    e_IV = sqrt(0.5)*randn(1,N);
    y_IV = lsim(B_true/A_true, u) + e_IV';
    
    x = lsim(B_LS/A_LS, u);
    
    % Regressor for instrumental variable
    % Zx
    Zx = zeros(N,n);    
    for iy = 1:N
        for jy = 1:n
            % Use at rest assumption
            if (iy-jy) < 1
                Zx(iy,jy) = 0;
            else
                Zx(iy,jy) = -x(iy-jy);
            end
        end
    end
    
    % Zu
    Zu = zeros(N,m);
    for iu = 1:N
        for ju = 1:m
            % Use at rest assumption
            if (iu-ju) < 1
                Zu(iu,ju) = 0;
            else
                Zu(iu,ju) = u(iu-ju);
            end
        end
    end
    
    Z = [Zx Zu];
    % Instrumental variable regression
    theta_IV = (Z'*Phi)\(Z'*y_IV);
    
    % A_IV
    A_IV = 1;
    for iy = 1:n
        A_IV = A_IV + theta_IV(iy) * z^(-iy);
    end
    
    % B_IV
    B_IV = 0;
    for iu = 1:m
        B_IV = B_IV + theta_IV(n+iu) * z^(-iu);
    end
    
    % Compare LS and IV method
    y_est_LS(:,k) = lsim(B_LS/A_LS, u);
    y_est_IV(:,k) = lsim(B_IV/A_IV, u);
    
    MSE_LS(k) = norm(y_true - y_est_LS(:,k))^2;
    MSE_IV(k) = norm(y_true - y_est_IV(:,k))^2;
    
    % Capture MSE which is too big
    if MSE_LS(k) > threshold
        MSE_LS(k) = inf;
    end
    if MSE_IV(k) > threshold
        MSE_IV(k) = inf;
    end 
end

% Plot MSE
figure(1);
plot(1:n_models, MSE_LS, 'ro--', 'MarkerSize', 10);
hold on
plot(1:n_models, MSE_IV, 'bo--', 'MarkerSize', 10);
title('MSE of LS and IV for different models');
legend({'LS', 'IV'});
xlabel('Model');
ylabel('MSE');
grid on

% Plot output
figure(2);
% Model 1
subplot(2,3,1);
plot(t,y_true,'r');
hold on
plot(t,y_est_LS(:,1),'b');
plot(t,y_est_IV(:,1),'g');
title('Model 1');
legend({'True', 'LS', 'IV'});
% Model 2
subplot(2,3,2);
plot(t,y_true,'r');
hold on
plot(t,y_est_LS(:,2),'b');
plot(t,y_est_IV(:,2),'g');
title('Model 2');
legend({'True', 'LS', 'IV'});
% Model 3
subplot(2,3,3);
plot(t,y_true,'r');
hold on
plot(t,y_est_LS(:,3),'b');
plot(t,y_est_IV(:,3),'g');
title('Model 3');
legend({'True', 'LS', 'IV'});
% Model 4
subplot(2,3,4);
plot(t,y_true,'r');
hold on
plot(t,y_est_LS(:,4),'b');
plot(t,y_est_IV(:,4),'g');
title('Model 4');
legend({'True', 'LS', 'IV'});
% Model 5
subplot(2,3,5);
plot(t,y_true,'r');
hold on
plot(t,y_est_LS(:,5),'b');
plot(t,y_est_IV(:,5),'g');
title('Model 5');
legend({'True', 'LS', 'IV'});

% Model 4 and 5 do not capture dynamics of true system and become unstable
% when using IV. IV is generally better then LS since this system is an
% output error model and IV leads to a smaller bias.