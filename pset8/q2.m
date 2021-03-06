% Setting parameters:

%parameters given
beta=0.99;
alpha=0.36;
delta=0.025;
k2y=10; % SS capital to output ratio
i2y=0.26; % SS investment to output ratio
rho=0.95;
sigma_e=0.007; % std dev of sigma
c2y=1-delta*k2y;

%steady state parameters
hBar=1/3;
A=(1-alpha)*(1-hBar)/(hBar*c2y);

rBar= 1/beta - (1-delta);
kBar=hBar*(alpha/rBar)^(1/(1-alpha));
yBar=kBar^alpha*hBar^(1-alpha);
cBar=yBar-delta*kBar;
wBar=(1-alpha)*yBar/hBar;
iBar=delta*kBar;
zBar=0;
B=wBar/cBar;

VARNAMES=['capital       '
          'consumption   ' 
          'output        ' 
          'investment    ' 
          'labor time    ' 
          'wages         ' 
          'rental rate   ' 
          'euler         ' 
          'technology    '];
      
AA=[0; 0; kBar; 0; 0; 0];
BB=[0; rBar*kBar; -kBar*(1-delta); alpha; alpha-1;alpha];
CC=[B*cBar 0 0 0 -wBar 0;
    -cBar 0 -iBar wBar*hBar wBar*hBar rBar*kBar;
    0 0 -iBar 0 0 0;
    0 -1 0 1-alpha 0 0;
    0 0 0 1-alpha 0 -1;
    0 0 0 -alpha -1 0 ];
DD= [0; 0; 0; 1; 1; 1];
FF=[0];
GG=[0];
HH=[0];
KK=[beta*(rBar+1-delta) 0 0 0 0 0];
JJ=[-beta*(rBar+1-delta) 0 0 0 0 beta*rBar];
NN=[rho];
MM=[0];
LL=[0];
Sigma=[sigma_e];

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);
HORIZON     = 32; % how far out should the impulse responses be calculated
PERIOD      = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
DO_PLOTS    = 0;  % if impulse response plots should be made, = 0, if not.
IMP_SELECT  = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
HP_SELECT   = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
GNP_INDEX   = 2; % Index of output among the variables selected for HP filter

do_it;

period=2000;
start=501;
time = 1:period;
error = sigma_e*randn(1,period);
z_t=zeros(1,period);
k_t=zeros(1,period);
Y_t=zeros(6,period);


for a = 2: period
    z_t(a)= rho*z_t(a-1)+error(a);
    k_t(a)=PP*k_t(a-1)+QQ*z_t(a);
    Y_t(:,a)=RR*k_t(a-1)+SS*z_t(a);
end
c_t=Y_t(1,start:period);
y_t_i=Y_t(2,start:period);
i_t=Y_t(3,start:period);
h_t_i=Y_t(4,start:period);
w_t=Y_t(5,start:period);
r_t=Y_t(6,start:period);

std_y=std(y_t);
vol_c=std(c_t)/std_y;
vol_i=std(i_t)/std_y;
vol_h=std(h_t)/std_y;
vol_w=std(w_t)/std_y;
vol_hw=std(h_t)/std(w_t);
corr_hw=corr(h_t',w_t');
% w=1600;
% [hpy,desvabsy] = hpfilter(y_t,w);
% [hph,desvabsy16] = hpfilter(h_t,w);
% cy=y_t-hpy';
% ch=h_t-hph';

figure(2)
scatter(y_t_i*100,h_t_i*100)
xlabel('output pct deviation from trend')
ylabel('hours worked pct deviation from trend')
title('Question 2:Indivisable Labor')