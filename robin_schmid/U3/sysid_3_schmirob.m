% System Identification Ex 3
% Robin Schmid, schmirob@ethz.ch
%% 2a
N = 1024;
e = randn(N,1); % Generate a random sequence (white Gaussian noise)
plot(e)
%% 2b
E = fft(e); % Fourier transform
%e_periodo = 1/N*abs(E).^2; % Do not use periodogram function, use own creation
%omega = linspace(0, 2*pi, N); % Frequencies from 0 to 2*pi

% Compare own creation to built in function,
% compare plots for omega from 0 to pi
e_periodo_2 = 1/N*abs(E(1:N/2+1)).^2;
e_periodo_builtin = periodogram(e);
omega_2 = linspace(0, pi, N/2+1);

loglog(omega_2, e_periodo_2, 'b');
hold on
loglog(omega_2, e_periodo_builtin, 'r');
factor = e_periodo_2./e_periodo_builtin
% Built in peridogram fct corresponds to 1/pi*e_periodo_2
%% 2c
N = 4096;
e = randn(N,1);
omega_2 = linspace(0, pi, N/2+1);

z = tf('z', 1); % T = 1 is the discrete sampling time
P = (z-0.6)/(z^2-0.4*z+0.85);
%P = tf([1, -0.6], [1, -0.4, 0.85], 1); % Define discrete plant directly

w = lsim(P, e); % Calculate linear system response
W = fft(w);
w_periodo = 1/N*abs(W(1:N/2+1)).^2;

% P for all frequencies, 1j: complex number
P = ((exp(1j.*omega_2)-0.6)./(exp(1j.*omega_2).^2-0.4.*exp(1j.*omega_2)+0.85))';

loglog(omega_2, abs(w_periodo(1:N/2+1)), 'g.-'); % Plot magnitude of output
hold on
loglog(omega_2, abs(P(1:N/2+1)), 'b.-'); % Plot magnitude of plant
loglog(omega_2, abs(w_periodo(1:N/2+1)-abs(P(1:N/2+1))), 'r.-'); % Plot magnitude of difference
xlabel('Frequency');
ylabel('Magnitude');
legend('Periodogram w','Magnitude P','Difference');