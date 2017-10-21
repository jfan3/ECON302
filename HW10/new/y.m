function income = y( j,k,t )

global z e W lambda
if t <= W
    income = exp(-2.17 + 0.168*(t+18) - 0.032*(t+18)^2/10 + 0.002*(t+18)^3/100 + z(j) + e(k));
else
    income = lambda*exp((-2.17 + 0.168*(W+18) - 0.032*(W+18)^2/10 + 0.002*(W+18)^3/100)+z(j));
end
