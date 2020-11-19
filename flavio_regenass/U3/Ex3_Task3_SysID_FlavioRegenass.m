%-------------------------------------------------
% System Identification Exercise Set 3, Task 3
% Flavio Regenass
% October 2020
%-------------------------------------------------
% Periodograms

%% Task 3, Part A
%----------------------------------------------------------------------------------------------------------------

% See the function file 'pseudoFFT.m'

%% Task 3, Part B
%% Subtask a)
% Generate e(k), 1024 point normally distributed random sequence
e = normrnd(0, 1, 4096,1);

%% Subtask b)
% Calculate Periodogram, PSD Estimate
[pxx, w] = periodogram(e);

% Check with Pseudo FFT
[E_N_self, w_self] = pseudoFFT(e);
pxx_self = 1/length(e) * abs(E_N_self).^2;

figure
plot(w, log10(pxx));
xlabel('Frequencies omega');
ylabel('log PSD for e_k');

%% Subtask c)
% Given Plant P(z) = (z - 0.6) / (z^2 - 0.4 * z + 0.85)
% Calculate the response w(k) of P to the input e(k)
syms z
P = (z - 0.6) / (z^2 - 0.4 * z + 0.85);

% Z Transform of Input Signal
y = sym('z');
e_Z = 0;
for i = 1:length(e)
    e_Z = e_Z + e(i)*y^(1-i);
end

% Output Signal in  Z-Domain
w_Z = P * e_Z;

% Inverse Transformation from Z-Domain to Discrete Domain to calculate w_k
sym_w_k = iztrans(w_Z);
w_k = zeros(length(e), 1);

for i = 1: length(e)
    n = i;
    w_k(i) = subs(sym_w_k);
end

%% Compute Omega_n for Plot of Transfer Function
omega_n = zeros(length(e), 1);
P_funcVal = zeros(length(e), 1);

for i = 1:length(e)
    omega_n(i) = 2*pi*(i-1)/length(e);
    P_funcVal(i) = (exp(sqrt(-1)*omega_n(i)) - 0.6) / (exp(sqrt(-1) * omega_n(i))^2 - 0.4 * exp(sqrt(-1) * omega_n(i)) + 0.85);
end

%% Subtask d)
% Periodogram of w_k

[wxx, v] = periodogram(w_k);
P_funcVal_trunk = abs(P_funcVal(1:((length(e)/2+1))));
error = abs(wxx) - P_funcVal_trunk;
error_squared = error.^2;

figure
loglog(v, wxx, v, abs(P_funcVal(1:(length(e)/2+1))).^2);
xlabel('Frequencies omega');
ylabel('Magnitude');
legend('Response of System W_N', 'Transferfunction Magnitude P', 'Error between the two')

figure
plot(v, sqrt(error_squared));
xlabel('Frequencies omega');
ylabel('Magnitude of Error');
legend('Error between Transfer Function P and Periodogram')
title('N = 4096')
