% Problem 4

load('experiment3.mat');

T = Camera_measurements.Temperature;
wavelength = Camera_sensitivity.Wavelenght; 

n = length(wavelength);
m = length(T);

emission_matrix = BlackBody_emission{:,2:end};
sens_R = Camera_sensitivity{:,2};
sens_G = Camera_sensitivity{:,3};
sens_B = Camera_sensitivity{:,4};

X_R = [];
X_G = [];
X_B = [];
for i=1:n
    X_R = [X_R sens_R(i)*emission_matrix(:,i)];
    X_G = [X_G sens_G(i)*emission_matrix(:,i)];
    X_B = [X_B sens_B(i)*emission_matrix(:,i)];
end

y_R = Camera_measurements{:,2};
y_G = Camera_measurements{:,3};
y_B = Camera_measurements{:,4};

% t_R = pinv(X_R)*y_R;

% close all
% plot(wavelength, t_R, 'x-')
% grid on

A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(n,1);
ub = ones(n,1);

t_R_constr = lsqlin(X_R,y_R,A,b,Aeq,beq,lb,ub);
t_G_constr = lsqlin(X_G,y_G,A,b,Aeq,beq,lb,ub);
t_B_constr = lsqlin(X_B,y_B,A,b,Aeq,beq,lb,ub);

figure(1)
close all
plot(wavelength, t_R_constr, 'x-', 'Color', 'r')
hold on
plot(wavelength, t_G_constr, 'x-', 'Color', 'g')
plot(wavelength, t_B_constr, 'x-', 'Color', 'b')
grid on
legend('Red', 'Green', 'Blue')
xlabel('\lambda')
ylabel('transmittance')


figure(2)
t_constr = lsqlin([X_R;X_G;X_B],[y_R;y_G;y_B],A,b,Aeq,beq,lb,ub);
plot(wavelength, t_constr, 'x-', 'Color', 'k')
grid on
xlabel('\lambda')
ylabel('transmittance')





