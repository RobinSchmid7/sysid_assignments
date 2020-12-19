% System Identification Ex 7
% Robin Schmid, schmirob@ethz.ch
%% 1.
clear all; close all; clc;

% Generate data
N_vec = 4.^(2:6); % N = [16 64 256 1024 4096]
N_exp = 1000;

theta1 = zeros(N_exp, 2);
theta2 = zeros(N_exp, 2);
it = 0;

for N = N_vec
    it = it + 1;
    
    f1 = figure(1);
    set(f1, 'visible', 'on');
    
    % Perform multiple experiments for statistical statement
    for e = 1:N_exp
        
        u = randn(N,1); % Gaussian with var = 1

        % a) Gaussian with var = 0.2
        w1 = sqrt(0.2)*randn(N,1);
        y1 = zeros(N,1);
        Phi1 = zeros(N, 2);

        y1(1) = w1(1);
        for i = 1:N-1
            y1(i+1) = 0.5*y1(i) + u(i) + w1(i+1);
        end

%         % Option 1: Assume system is initially at rest, use full y vector with
%         % data
%         Phi1(1,:) = [0 0];
%         for k = 2:N
%             Phi1(k,:) = [y1(k-1) u(k-1)];
%         end
% 
%         % Option 2: Assume data at t = -1 is the same as on t = 0
%         Phi1(1,:) = [y1(1) u(1)];
%         for k = 2:N
%             Phi1(k,:) = [y1(k-1) u(k-1)];
%         end

        % Option 3: Discard negative data, here discard y(1)
        % Option 3 and 1 lead to the same solution
        for k = 1:N-1
            Phi1(k,:) = [y1(k) u(k)];
        end
        Phi1(N,:) = [];
        y1(1) = [];

        % LS solution
        theta1(e,:) = Phi1\y1;

        % b) Uniform with var = 0.2
        b = sqrt(6*0.2); % Half of size of interval
        w2 = 2*b*rand(N,1) - b;
        y2 = zeros(N,1);
        y2(1) = w2(1);
        for i = 1:N-1
            y2(i+1) = 0.5*y2(i) + u(i) + w2(i+1);
        end
        
        % Option 1: Assume system is initially at rest, use full y vector with
        % data
        Phi2(1,:) = [0 0];
        for k = 2:N
            Phi2(k,:) = [y2(k-1) u(k-1)];
        end 
        
        % LS solution
        theta2(e,:) = Phi2\y2;
        
    end
    
    % Plotting
    legend_info{it} = num2str(N);
    
    subplot(2,2,1);
    histogram(theta1(:,1));
    title('Histogram for parameter a, Gaussian noise');
    xlabel('a');
    ylabel('Number of samples');
    legend(legend_info);
    hold on
    grid on

    subplot(2,2,3);
    histogram(theta1(:,2));
    title('Histogram for parameter b, Gaussian noise');
    xlabel('b');
    ylabel('Number of samples');
    legend(legend_info);
    hold on
    grid on
    
    subplot(2,2,2);
    histogram(theta2(:,1));
    title('Histogram for parameter a, uniform noise');
    xlabel('a');
    ylabel('Number of samples');
    legend(legend_info);
    hold on
    grid on

    subplot(2,2,4);
    histogram(theta2(:,2));
    title('Histogram for parameter b, uniform noise');
    xlabel('b');
    ylabel('Number of samples');
    legend(legend_info);
    hold on
    grid on
    
    % For bigger N the variance of the estimate decreases, effect of noise
    % gets filtered out, method a) is better than b), i.e. Gaussian noise
    % is better than uniform noise, i.e. effect is reduced more when using
    % a LS estimate
end

%% 2.
% Generate data
N_vec = 4.^(3:6); % N = [64 256 1024 4096]
N_exp = 1000;

theta1 = zeros(N_exp, 2);
it = 0;
total_err = zeros(N_exp,1);

f2 = figure(2);
set(f2, 'visible', 'on');

for N = N_vec
    it = it + 1;
    
    % Perform multiple experiments for statistical statement
    for e = 1:N_exp
        
        u = randn(N,1); % Gaussian with var = 1

        % a) Gaussian with var = 0.2
        w1 = sqrt(0.2)*randn(N,1);
        y1 = zeros(N,1);
        Phi1 = zeros(N, 2);

        y1(1) = w1(1);
        for i = 1:N-1
            y1(i+1) = 0.5*y1(i) + u(i) + w1(i+1);
        end

%         % Option 1: Assume system is initially at rest, use full y vector with
%         % data
%         Phi1(1,:) = [0 0];
%         for k = 2:N
%             Phi1(k,:) = [y1(k-1) u(k-1)];
%         end

%         % Option 2: Assume data at t = -1 is the same as on t = 0
%         Phi1(1,:) = [y(1) u(1)];
%         for k = 2:N
%             Phi1(k,:) = [y1(k-1) u(k-1)];
%         end

        % Option 3: Discard negative data, here discard y(1)
        for k = 1:N-1
            Phi1(k,:) = [y1(k) u(k)];
        end
        Phi1(N,:) = [];
        y1(1) = [];

        % LS solution
        theta1(e,:) = Phi1\y1;

        % Sum of squared residual
        total_err(e) = norm(y1 - Phi1 * theta1(e,:)')^2;
    end
    
    % Normalizing by variance
    total_err(:) = total_err(:) / 0.2;
    
    % Plotting
    legend_info{it} = num2str(N);
    
    subplot(2,2,it);
    histogram(total_err);
    hold on
    grid on
    histogram(N-2 + sqrt(2*(N-2))*randn(N_exp,1)); % Gaussian
    histogram(chi2rnd(N-2,[N_exp,1])); % Chi square
    xlabel('|e|^2');
    ylabel('Number of samples');
    legend({legend_info{it}, 'normal Gaussian', 'chi square'});
    
    % Chi square fits to sum of squared residual distribution, Gaussian
    % distribution is a good approximation of chi squared for big N
end