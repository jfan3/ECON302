
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
k_newStar= k_star-0.1*g_L;
k_vect = 0:0.1:35;
c_lowTax = k_vect.^alpha-g_L-k_vect.*(n+gama+delta);
c_highTax= k_vect.^alpha-g_H-k_vect.*(n+gama+delta);
c_lowStar = k_star^alpha-g_L-k_star*(n+gama+delta);
c_highStar= k_newStar^alpha-g_H-k_newStar*(n+gama+delta);

figure (1)
plot(k_vect,c_lowTax,'b',k_vect,c_highTax,'r',[k_star k_star],[-0.6 1.2],'k')
hold on
scatter([k_star k_star],[c_lowStar c_highStar],'r')