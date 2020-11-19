% This function takes the vector u and transforms it into the FFT-like
% Vector U_N

function [ U_N, omega ] = pseudoFFT(u)

% check if only real values
check = isreal(u);

N = length(u);

omega = zeros(N, 1);
U_N = zeros(N, 1);

% Adjust Vector length if Real input vector
if check == 1
    %N = ceil((N+1)/2);
    
    omega = zeros(N, 1);
    U_N = zeros(N, 1);
    
    for k = 1:N
        % Fill FFT Array U_N(k)
        for i = 1:N
            omega(i) = 2*pi*(i-1)/N; % Caution here with the indexing, since Matlab start at 0, but omega starts with omega_0, where n = 0 = i
            U_N(k) = U_N(k) +  u(i) * exp(-sqrt(-1) * omega(i)*(k-1)); % Notice the factor two, since the FFT is only composed of complex conjugate pairs
        end
    end
    
elseif check == 0
    
    omega = zeros(N, 1);
    U_N = zeros(N, 1);
    
    for k = 1:N
        % Fill FFT Array U_N(k)
        for i = 1:N
            omega(i) = 2*pi*(i-1)/N; % Caution here with the indexing, since Matlab start at 0, but omega starts with omega_0, where n = 0 = i
            U_N(k) = U_N(k) + u(i) * exp(-sqrt(-1) * omega(i)*(k-1));
        end
    end
end

% realU_N = fft(u);
% 
% error = U_N-realU_N;


end

