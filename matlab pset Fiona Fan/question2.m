%%Question 2
% The function for question a is defined as function "qu_a"

%Definition
K = 0.1:0.1:10;
index = 1;
%placeholders
diff = zeros (1,length(K));
pig = zeros (1,length(K));

%diff is the the difference between the marginal 
%product of capital and the rental rate of capital
%pig is the profit
while index<=length(K)
diff(index)=qu_a(K(index));
pig(index) = K(index)^(1/3)*10^(2/3)-K(index)-2;
index=index+1;
end

%Question 2b: Plot the difference and profit
figure (1)
plot(K, diff,'g')
hold on
plot(K, pig, 'r')
xlabel ('K($)')
ylabel ('$')
legend('MPC/price difference','profit')

%negative profit
neg_prof = @(K)(K^(1/3)*10^(2/3)-K-2)*(-1);
%To find the point where the marginal product of capital equals the rental rate of capital
FOC_0=fzero(@qu_a,10) %1.9245
FOC_1=bisect(@qu_a,0,10) %using the bisection method,1.9238
%Find the minmum of the negatie profit
min = fminunc(neg_prof,3) %1.9245