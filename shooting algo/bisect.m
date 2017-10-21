% Bisection method to find zero. Need the function, guestimate of lower and
% upper bound that encompass the zero point of the function. Only suitable
% for one zero point
function x0 = bisect (func, lower, upper)
%if the lower and upper bounds are on the same size of x axis, then there's
%no point crossing zero within the bound.
if func(lower)*func(upper)>0
    x0='N/A';
else 
    mid = (lower+upper)/2;
    while abs(func (mid))>10^(-3)
        %If the mid point is between lower bound the mid point, update
        %upper bound to mid point. Vice versa.
        if func(lower)*func(mid)<0
            upper = mid;
        else
            lower = mid;
        end
    mid = (lower+upper)/2;
    
    end
    %If the value of mid point is smalelr than 1*10^-3, consider the mid
    %point as the zero point
    x0 = mid;
end