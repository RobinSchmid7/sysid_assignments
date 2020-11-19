function [R_rand,lags,phi_u,omegas] = HS2020_SysID_Exercise_05_16927725
%% Output format specification
% R_u must be a 2xN matrix
% phi_u must be a 2xN matrix
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[u_prbs,u_rand] = HS2020_SysID_Exercise_05_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the input signals u_prbs and u_randn to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% 1. Calculation of autocorrelation
N = length(u_prbs); % Here calculation length = period
% Option 1: Compute only positive lags and mirror
%lags = [0:(N-1)/2];
% Option 2: Compute all lags
lags = [-(N-1)/2:(N-1)/2];

% Option 1: Calculate R_prbs exactly
% u = u_prbs(1).^2;
% R_prbs(1) = u;
% R_prbs(2:N) = -u/N;
% plot(lags, R_prbs, 'b',  -lags, R_prbs, 'b');

% Option 2: Via definition of autocorrelation
n_lags = size(lags, 2);
R_prbs = zeros(1, n_lags);
for t = 1:n_lags
    for k = 1:N
        % Account for periodicity of signal
        if k-lags(t)<=0
            R_prbs(t) = R_prbs(t) + 1/N * u_prbs(k) * u_prbs(k-lags(t)+N);
        elseif k-lags(t)>N
            R_prbs(t) = R_prbs(t) + 1/N * u_prbs(k) * u_prbs(k-lags(t)-N);
        else
            R_prbs(t) = R_prbs(t) + 1/N * u_prbs(k) * u_prbs(k-lags(t));
        end
    end
end

R_rand = zeros(1, n_lags);
% Would also need to account for signals with an uneven number of samples,
% lags must be integer numbers, in this case calculate from 0 to N
for t = 1:n_lags
    for k = 1:N
        % Account for periodicity of signal
        if k-lags(t)<=0
            R_rand(t) = R_rand(t) + 1/N * u_rand(k) * u_rand(k-lags(t)+N);
        elseif k-lags(t)>N
            R_rand(t) = R_rand(t) + 1/N * u_rand(k) * u_rand(k-lags(t)-N);
        else
            R_rand(t) = R_rand(t) + 1/N * u_rand(k) * u_rand(k-lags(t));
        end
    end
end

% Plotting
f1 = figure(1);
set(f1, 'visible', 'off');
p1 = plot(lags, R_prbs, 'b');
hold on
% plot(-lags, R_prbs, 'b'); % For mirroring with only positive lags
p2 = plot(lags, R_rand, 'r');
% plot(-lags, R_rand, 'r'); % For mirroring with only positive lags
legend([p1 p2], {'R_{prbs}', 'R_{rand}'});
xlabel('lag');
ylabel('R_u');

R_u = [R_prbs; R_rand];
%% Calculate spectra
omegas = (2*pi/N)*[0:N-1]';

% Option 1: Calculate phi_prbs exactly
% u = u_prbs(1).^2;
% phi_prbs(1) = u*(1-(N-1)/N);
% phi_prbs(2:N) = u*(1+1/N);

% Option 2: Calculate via fourier transform
phi_prbs = fft(R_prbs);
phi_rand = fft(R_rand);

% Plotting
f2 = figure(2);
set(f2, 'visible', 'on');
semilogx(omegas, abs(phi_prbs), 'b'); % Plot magnitude
hold on
semilogx(omegas, abs(phi_rand), 'r'); % Plot magnitude
xlim([1e-2,max(omegas(N))]);
ylim([0,1.2]);
legend({'phi_{prbs}', 'phi_{rand}'});
xlabel('lag');
ylabel('phi_u');

phi_u = [phi_prbs; phi_rand];
end