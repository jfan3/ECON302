%%Question 1

%x is a column vector from 0.01 to 1
x = (0.01:0.01:1)';
unos = ones (100,1);

%epsilon is randomly drawn from normal distribution
epsilon = randn (100,1);

%y = 10x + 5 + epsilon
y = [unos x]*[5;10]+epsilon;
%X = [1 x]
x_big = [unos x];

% beta estimated = (X'*X)^(-1)*X'*Y
% intercept ~= 5, slope~=10
beta_e = inv((x_big')*x_big)*x_big'*y

%calculate estimated y
y_estimated =x_big*beta_e;

figure(1)
subplot(2,1,1)
plot (x,y_estimated,'r')
hold on
plot (x,y)
title ('OLS estimation')
xlabel('x')
ylabel('y')
legend('estimated y', 'y')


residuals = y-x_big*beta_e;
subplot(2,1,2)
hist (residuals,25)
title ('residual distribution')
xlabel('residual')
ylabel('frequency')





