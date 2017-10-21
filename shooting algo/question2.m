%global variables that never change
global T beta gama n delta alpha A0 N0;
T =100;
beta=0.96;
gama=0.02;
n=0.01;
delta=0.06;
alpha=0.36;
A0=1;
N0=1;

%initial conditions
k_star= (alpha/(beta^(-1)*(1+gama)-(1-delta)))^(1/(1-alpha));
y_star= k_star^alpha;
g_L= 0.2*y_star;
g_H= 1.1*g_L;
k_newStar = k_star-0.1*g_L;
% k0=k_star;
k_series= zeros(1,T);
k_series(1)=k_star;
lower = 0;
upper = k_star;


i=1;
difference = k_star-0;

t_series = ones(1,100)*g_L;
t_series(1:5)=ones(1,5)*g_H;

k_vect = 0:0.1:35;
c_lowTax = k_vect.^alpha-g_L-k_vect.*(n+gama+delta);
c_highTax= k_vect.^alpha-g_H-k_vect.*(n+gama+delta);
c_lowStar = k_star^alpha-g_L-k_star*(n+gama+delta);
c_highStar= k_newStar^alpha-g_H-k_newStar*(n+gama+delta);

ki = (lower+upper)/2;
%question a. Shooting algorithm in its flesh
while abs(difference)>0.001
    %While kt hasn't arrived at k*, guess a k1, use k1 and k0 as the
    %initial conditions for second order differentiation to propagate a
    %series of k. 
    k_series(3:T)= zeros(1,T-2);
    
    k_series(2)=ki;
    
    %For one k1, when k_t reaches towards 0 or bigger than k_star, exit the
    %propagation and start with a new k1.
    for a=3:T
    k_series(a)=sdeq(k_series(a-2),k_series(a-1),t_series(a-2),t_series(a-1),alpha, beta,delta,gama,n);
        if k_series(a)<0 ||k_series(a)>k_star
            break
        end
    end
    
    %Use bisection method to guess k1. Update the upper and lower bounds of
    %the guess.
    if k_series(a)>k_star
        upper=ki;
            
    else 
        lower=ki;
    end
    
    i=i+1;
    %Update the difference between the last element in the series and k*
    difference = k_series(T)-k_star;
    ki = (lower+upper)/2;
end

%question b-f
c_series= zeros(1,T-1);
kStagnant=zeros(1,T);
fkt_series=k_series.^alpha;
cStagnant=ones(1,T).*k_star;
r_series = zeros(1,T);
i_series = zeros(1,T-1);
sr_series= zeros(1,T-1);
for b = 1:T-1
    c_series(b)=k_series(b)*(1-delta)+k_series(b)^alpha-k_series(1+b)*(1+gama)*(1+n)-t_series(b);
    i_series(b) = fkt_series(b)-c_series(b);
    sr_series(b) = i_series(b)/fkt_series(b);
end
for c = 1:T
    kStagnant(c)=fkt_series(c)-(n+gama+delta)*k_series(c);
    r_series(c) = alpha*k_series(c)^(alpha-1);
    
end

time1 = 1:T-1;
time2 = 1:T;

A_series=zeros(1,T);
N_series=zeros(1,T);
Y_series=zeros(1,T);
C_series=zeros(1,T-1);
K_series=zeros(1,T);
I_series=zeros(1,T-1);
for d = 1:T
    A_series(d)=A0*(1+gama)^(d-1);
    N_series(d)=N0*(1+n)^(d-1);
    Y_series(d)=fkt_series(d)*A_series(d)*N_series(d);
    K_series(d)=k_series(d)*A_series(d)*N_series(d);
end
for e =1:T-1
    C_series(e)=c_series(e)*A_series(e)*N_series(e);
    I_series(e)=i_series(e)*A_series(e)*N_series(e);
end

figure (1)
plot(k_vect,c_lowTax,'b',k_vect,c_highTax,'r',[k_star k_star],[-0.6 1.2],'k')
legend('\Deltak=0_b','\Deltak=0_a', '\Deltac=0')

figure (2)
plot(k_vect,c_lowTax,'b',k_vect,c_highTax,'r',[k_star k_star],[-0.6 1.2],'k')
% legend('c_before','c_after', '\Deltac=0')
hold on
scatter(k_series(1:T-1),c_series)
legend('\Deltak=0_b','\Deltak=0_a', '\Deltac=0','Evolution')
axis([4.5 6.5 0.9 1])