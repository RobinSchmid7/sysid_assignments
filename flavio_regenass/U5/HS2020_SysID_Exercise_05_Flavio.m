function [R_u,lags,phi_u,omegas] = HS2020_SysID_Exercise_05_16914186
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

%% 1. Calculation of autocorrelation ( Slide 4.20 for generic Signal and Slide 6.48 for PRBS)
Ts = 1; % Sampling Time
N = length(u_prbs); % Calculation Length = Period here

R_u = zeros(2, N);
lags = linspace(-(N-1)/2, (N-1)/2, 127); % Lags for Autocorrelation

for tau = -(N-1)/2:1:(N-1)/2
    sum_autocorr = [0, 0]; % 1x2 Vector for the two Input signals u_prbs and u_rand
    for k = 1:N
        if tau == 0
            sum_autocorr = sum_autocorr + [u_prbs(k)^2, u_rand(k)^2];
        elseif (k-tau >= 1 && k-tau <= N)
            sum_autocorr = sum_autocorr + [u_prbs(k)*u_prbs(k-tau), u_rand(k)*u_rand(k-tau)];
        elseif (k-tau < 1) % Since we assume our Signals to be periodic before and after the measurement u(k) = u(k + N)
            sum_autocorr = sum_autocorr + [u_prbs(k)*u_prbs(k+N-tau), u_rand(k)*u_rand(k+N-tau)];
        elseif (k-tau > N)
            sum_autocorr = sum_autocorr + [u_prbs(k)*u_prbs(k-tau-N), u_rand(k)*u_rand(k-tau-N)];
        end
    end
    R_u(:, (N-1)/2+1 + tau) =  1/N .* sum_autocorr.';
    
%    Direct Approach for PRBS Signal
%     if tau == 0
%         R_u(1, (N-1)/2+1) = max(u_prbs)^2;
%     else
%         R_u(1, (N-1)/2+1 + tau) = -max(u_prbs)^2/N;
%     end
end


figure(1);
plot(lags, R_u(1, :), '--x' , lags, R_u(2,:), ':o');
legend('R_{u}_{prbs}', 'R_{u}_{rand}');
xlabel('Lags Tau');
title('Autocorrelation of Inputs u_{prbs} and u_{rand}');

figure(2);
plot(linspace(1,N,127), u_prbs, linspace(1, N, 127), u_rand);
legend('Input PRBS', ' Input Rand');
title('Inputs u_{prbs} and u_{rand}');

% R_u = NaN*ones(2,N);
% lags = NaN*ones(1,N);
%% Calculate spectra ( Slide 4.21 for generic Signal, 6.50 for PRBS)

% First Attempt:
omega = 2*pi/N * [0:N-1]; % Frequencies for PSD
% phi_u = zeros(2,N);
% % for n = 1:N
% %     for tau = (N-1)/2:N
% %         phi_u(:, tau) = phi_u(:, tau) + R_u(:, tau) .* exp(-(sqrt(-1))*omega(n)* (tau-1));
% %     end
% % end
% 
phi_fft = fft(R_u, [], 2);
phi_fft_prbs = fft(R_u(1,:), 127);
phi_fft_rand = fft(R_u(2,:), 127);
delta_1 = phi_fft_prbs - phi_fft(1,:)
delta_2 = phi_fft_rand - phi_fft(2,:)
phi_spectr_fft = phi_fft;
% [phi_PRBS_periodogram, w_1] = periodogram(u_prbs);
% [phi_RAND_periodogram, w_2] = periodogram(u_rand);
% 
% figure(3);
% plot(w_1, phi_PRBS_periodogram, w_2, phi_RAND_periodogram);
% legend('Periodogram PRBS', 'Periodogram Rand');
% figure(4);
plot(omega(1:(N-1)/2+1), abs(phi_fft_prbs(1:(N-1)/2+1)), omega(1:(N-1)/2+1), abs(phi_fft_rand(1:(N-1)/2+1)));
legend('Spectrum using FFT(Autocorrelation(PRBS))', 'Spectrum using FFT(Autocorrelation(Rand))');
title('FFT Approach for Spectrum Calculation');
% 
% figure(5);
% plot(omega, phi_fft(1,:), omega, phi_fft(2, :));
% legend('FFT(PRBS)', 'FFT(Rand)');
% 
% phi_u = [phi_PRBS_periodogram; ...
%         phi_RAND_periodogram];
% omegas = w_1

% Calculate FFT of Input Signals
% (https://ch.mathworks.com/help/matlab/math/fft-for-spectral-analysis.html)
U_prbs = fft(u_prbs, 127);
U_rand = fft(u_rand, 127);
% Calculate Spectrum from Absolute Value (Slide 4.18)
phi_u_prbs = 1/N * U_prbs.*conj(U_prbs)/127;
phi_u_rand = 1/N * U_rand.*conj(U_rand)/127;
% Frequencies of FFT, we only use half of the frequencies since the rest
% will be the complex conjugate
omegas = 2*pi/127*(0:126);
% Plot
figure(6);
plot(omegas(1:66), phi_u_prbs(1:66), omegas(1:66), phi_u_rand(1:66));
xlabel('Freq [rad/s]');
ylabel('Spectrum of inputs u');
legend('Phi from u_{prbs}', 'Phi from u_{rand}');
title('With Solution from Matlab Help Page');

phi_u = [phi_u_prbs; phi_u_rand];
omegas = omegas;
%phi_u = NaN*ones(2,N);
%omegas = NaN*ones(1,N);
end
