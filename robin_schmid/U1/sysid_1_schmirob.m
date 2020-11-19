%System Identification Ex 1
%Author: Robin Schmid, schmirob@ethz.ch
%% 1 (maximum likelihood)
w = u + 2*u.^2 - u.^3;
A = w'*w;
b = (y - 1.5 * ones(100, 1))'*w;
x = A\b;
%% 2 (maximum a posterior)
w = u + 2*u.^2 - u.^3;
A = 1/1.2 * (w'*w) + 1/0.09;
b = 1/1.2 * (y - 1.5 * ones(100,1))'*w + 1.5/0.09;
x = A\b
