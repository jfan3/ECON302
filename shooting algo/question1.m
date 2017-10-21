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
k0=0.5*k_star;

k_series= zeros(1,T);
k_series(1)=k0;
lower = k0;
upper = k_star;
% ki = (lower+upper)/2;
i=1;
difference = k_star-k0;
d_series=zeros(5000,1);
d_series(1)=difference;
g_L= 0.2*y_star;
g_H= 1.1*g_L;
k_newStar=k_star-0.1*g_L;


%question a. Shooting algorithm in its flesh
while abs(difference)>0.001
    %While kt hasn't arrived at k*, guess a k1, use k1 and k0 as the
    %initial conditions for second order differentiation to propagate a
    %series of k. 
    k_series(3:T)= zeros(1,T-2);
    ki = (lower+upper)/2;
    k_series(2)=ki;
    
    %For one k1, when k_t reaches towards 0 or bigger than k_star, exit the
    %propagation and start with a new k1.
    for a=3:T
    k_series(a)=sdeq(k_series(a-2),k_series(a-1),alpha, beta,delta,gama,n);
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
    d_series(i)=difference;
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
    c_series(b)=k_series(b)*(1-delta)+k_series(b)^alpha-k_series(1+b)*(1+gama)*(1+n);
    i_series(b) = fkt_series(b)-c_series(b);
    sr_series(b) = i_series(b)/fkt_series(b);
end
for c = 1:T
    kStagnant(c)=fkt_series(c)-(n+gama+delta)*k_series(c);
    r_series(c) = alpha*k_series(c)^(alpha-1);
    
end

time1 = 1:T-1;
time2 = 1:T;

figure(1)
title('Phase diagram')
plot(k_series(1:T-1),c_series,'b',k_series,kStagnant,'r',[k_star k_star],[0 2],'k')
legend('k_t vs c_t','\Deltak=0','\Deltac=0')
xlabel('k_t')
ylabel('c_t')

figure(2)
title('rates')
subplot(1,2,1)
plot(time2, r_series,'k',time1,sr_series,'m')
legend('rental rate','savings rate')
xlabel('time')
ylabel('rate')
subplot(1,2,2)
plot (sr_series, r_series(1:T-1))
xlabel('saving rate')
ylabel('rental rate')

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
figure(3)
title('Aggregates')
subplot(2,2,1)
plot(time2,Y_series)
xlabel('time')
ylabel('Y')
subplot(2,2,2)
plot(time1,C_series)
xlabel('time')
ylabel('C')
subplot(2,2,3)
plot(time2,K_series)
xlabel('time')
ylabel('K')
subplot(2,2,4)
plot(time1,I_series)
xlabel('time')
ylabel('I')

figure(4)
title('Exploring differences')
plot(1:i,d_series(1:i))
xlabel('i')
ylabel('difference')