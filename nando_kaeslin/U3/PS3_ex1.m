% System Identification - Exercise 3 MATLAB, Problem 1

close all

t = 0:0.01:pi;
u = zeros(length(t),1);

for i=20:10:60
    for j=1:length(u)
        u(j) = u(j)+sin(i*t(j));
    end
end

figure(1)
plot(t,u)
title('Signal in Time domain')
xlabel('t')
ylabel('u(t)')

grid on

N = length(u);
omega_n = 0:2*pi/N:2*pi*(N-1)/N;

U = fast_fourier_transform(u, omega_n);

figure(2)
plot(omega_n, abs(U))
grid on
title('Fourier Transform (my implemetation)')
xlabel('\omega')
ylabel('U(\omega)')

figure(3)
U = fft(u);
plot(abs(U))
grid on
title('Fourier Transform (fft( ) implemetation)')
xlabel('\omega')
ylabel('U(\omega)')


function U_N = fast_fourier_transform(u, omega_n)

    N = length(u);

    U_N = zeros(1,N);
    for i=1:N
        for j=1:N
            U_N(i) = U_N(i)+u(j)*exp(-1i*omega_n(i)*(j-1));
        end
    end
end



