%Autocorrelation exercise
%Author: Robin Schmid, schmirob@ethz.ch
%%
N = 400;
t = 0:N-1;
%u = 1.*rand(N,1);
%u = sawtooth(t.*2*pi/50, 1/2)';
u = idinput(N);
lags = [0:((N-1)/2)]; % For positive time lags
% figure (1)
% plot(t, u, 'b')

r = auto_corr(u, lags, N);
figure (2)
plot(lags,r,'b')
hold on
plot(-lags,r,'b') % Plot mirror around y axis

function [R_u] = auto_corr(u, lags, N)
    n_lags = size(lags, 2);
    R_u = zeros(1, n_lags);
    for t = 1:n_lags
        for k = 1:N
            % Account for periodicity of signal
            if k-lags(t)<=0
                R_u(t) = R_u(t) + 1/N * u(k) * u(k-lags(t)+N);
            elseif k-lags(t)>N
                R_u(t) = R_u(t) + 1/N * u(k) * u(k-lags(t)-N);
            else
                R_u(t) = R_u(t) + 1/N * u(k) * u(k-lags(t));
            end
        end
    end
end





