global T beta gama n delta alpha A0 N0;
T =1000;
beta=0.96;
gama=0.02;
n=0.01;
delta=0.06;
alpha=0.36;
A0=1;
N0=1;

k_star= (alpha/(beta^(-1)*(1+gama)-(1-delta)))^(1/(1-alpha));
k0=0.5*k_star;

k_series= zeros(1,T);
k_series(1)=k0;
lower = k0;
upper = k_star;
% ki = (lower+upper)/2;
i=1;

difference = k_star-k0;

while abs(difference)>0.001
   
     k_series(3:T)= zeros(1,T-2);
    ki = (lower+upper)/2;
    k_series(2)=ki;
    
    for a=3:100
    k_series(a)=sdeq(k_series(a-2),k_series(a-1),alpha, beta,delta,gama,n);
        if k_series(a)<0 ||k_series(a)>k_star
            break
        end
    end
    
    
    if k_series(a)>k_star
        upper=ki;
            
    else 
        lower=ki;
    end
    
    i=i+1;
    difference = k_series(T)-k_star;
end

c_series= zeros(1,T-1);
kStagnant=zeros(1,T);
fkt_series=k_series.^alpha;
cStagnant=ones(1,T).*k_star;
r_series = zeros(1,T);
i_series = zeros(1,T-1);
sr_series= zeros(1,T-1);
for b = 1:99
    c_series(b)=k_series(b)*(1-delta)+k_series(b)^alpha-k_series(1+b)*(1+gama)*(1+n);
    i_series(b) = fkt_series(b)-c_series(b);
    sr_series(b) = i_series(b)/fkt_series(b);
end
for c = 1:100
    kStagnant(c)=fkt_series(c)-(n+gama+delta)*k_series(c);
    r_series(c) = alpha*k_series(c)^(alpha-1);
    
end

time1 = 1:99;
time2 = 1:100;

figure(1)
% plot(time1, c_series,'b',time2,k_series,'g',time2,kStagnant,'r', time2,cStagnant,'k')
% legend('consumption','capital','\Deltak=0','\Deltac=0')
% xlabel('time')
% ylabel('output')
plot(k_series,c_series,'b',k_series,kStagnant,'r',[k_star k_star],[0 2],'k')
legend('k_t vs c_t','\Deltak=0','\Deltac=0')
xlabel('k_t')
ylabel('c_t')

figure(2)
subplot(1,2,1)
plot(time2, r_series,'k',time1,sr_series,'m')
legend('rental rate','savings rate')
xlabel('time')
ylabel('rate')
subplot(1,2,2)
plot (sr_series, r_series(1:99))
xlabel('saving rate')
ylabel('rental rate')

A_series=zeros(1,100);
N_series=zeros(1,100);
Y_series=zeros(1,100);
C_series=zeros(1,99);
K_series=zeros(1,100);
I_series=zeros(1,99);
for d = 1:100
    A_series(d)=A0*(1+gama)^(d-1);
    N_series(d)=N0*(1+n)^(d-1);
    Y_series(d)=fkt_series(d)*A_series(d)*N_series(d);
    K_series(d)=k_series(d)*A_series(d)*N_series(d);
end
for e =1:99
    C_series(e)=c_series(e)*A_series(e)*N_series(e);
    I_series(e)=i_series(e)*A_series(e)*N_series(e);
end
figure(4)
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

