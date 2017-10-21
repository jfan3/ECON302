% Quiz #2 Part 1
% 11/11/2016
% Fiona Fan && Adrian Cao

%(f)
clear all

data = xlsread('Data.xlsx');
[length num_series]=size (data);

Cn=data(:, 2);
Cd=data(:, 3);
i=data(:, 4);
% Average ratios
avgCnCd=mean(Cn./Cd);
avgCni=mean(Cn./i);
avgCdi=mean(Cd./i);
y=log(data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodrick-Prescott Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the trend component of the data, also in logs
y_HP = zeros(length,num_series);
for i = 1:num_series
    [y_HP(:,i)] = hpfilter(y(:,i),1600);
end

% Compute filtered series of residuals in percentage deviation from trend
u_HP = y - y_HP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Correlation Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CorrMatrix = corr(u_HP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Relative Volatilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RelVol = std(u_HP)/std(u_HP(:,1));

% Table
table=[RelVol; CorrMatrix];
% The ratio of nondurable goods consumption to durable goods consumption is 
% much 6.511, meaning that there are more nondurable goods in the market
% than durable ones.
% The consumption of durable goods is much more volatile than that of
% nondurable goods when compare to output (2.93 VS 0.54). This makes sense 
% because people are more elastic to durable goods consumption because they
% cost a lot, while for non-durable goods like food and clothes, they are
% less elastic. 

% (g)
% Given values
% Fiona's initial values
% z_bar     = 1; % Normalization
alpha     = 0.36; % Capital's share of total income
delta     = 0.025; % Depreciation rate for capital
rho      = 0.95; % autocorrelation of technology shock
sigma_e = 0.7; % Standard deviation of technology shock.  Units: Percent.
% 
beta     = 0.99;
A=(1-beta*(1-delta))/(delta*beta*avgCnCd);
RBar     = 1/beta-(1-delta); % steady state return on capital before depreciation
PkBar = 1;
PdBar=1;
WBar=(1-alpha)/alpha*((alpha*(RBar^(-alpha)))^(1-alpha));
hk_ratio= RBar*(1-alpha)/(WBar*alpha);
rw_ratio=alpha/(1-alpha)*hk_ratio;
HkdBar=0.18;
HddBar=0.11;
HndBar=0.71;
CnBar=hk_ratio^(-alpha)*HndBar;
CdBar=hk_ratio^(-alpha)*HddBar;
IBar=hk_ratio^(-alpha)*HkdBar;
KndBar=CnBar/(hk_ratio^(1-alpha));
KddBar=CdBar/(hk_ratio^(1-alpha));
KkdBar=IBar/(hk_ratio^(1-alpha));
KBar=KndBar+KddBar+KkdBar;
DBar=CdBar/delta/PdBar;
ZBar=0;
YBar=CnBar+CdBar+IBar;
% 
% %Adrian's initial values
% alpha=0.36;
% beta=0.99;
% delta=0.025;
% rho=0.95;
% sigma_e=0.7;
% % Calculate A
% A=(1-beta*(1-delta))/(6.5111*delta*beta);
% HkdBarKkdBar=((1-beta*(1-delta))/(alpha*beta))^(1/(1-alpha));
% % Steady State
% ZBar=0;
% KBar=1/HkdBarKkdBar;
% KkdBar=delta*KBar/(HkdBarKkdBar^(1-alpha));
% HkdBar=HkdBarKkdBar*KkdBar;
% DBar=avgCdi*KBar;
% KddBar=delta*DBar/HkdBarKkdBar;
% KndBar=KBar-KkdBar-KddBar;
% HddBar=HkdBarKkdBar*KddBar;
% HndBar=HkdBarKkdBar*KndBar;
% RBar=alpha*KndBar^(alpha-1)*HndBar^(1-alpha);
% WBar=(1-alpha)*KndBar^(alpha-1)*HndBar^(1-alpha);
% PkBar=beta*RBar/(1-beta*(1-delta));
% PdBar=PkBar;
% CdBar=delta*DBar*PdBar;
% CnBar=KndBar^alpha*HndBar^(1-alpha);
% IBar=delta*KBar*PkBar;
% YBar=CnBar+CdBar+IBar;

% (h)
VARNAMES=['k ',
          'd ',
          'cn',
          'cd',
          'i ',
          'y ',
          'r ',
          'w ',
          'kn',
          'kd',
          'kk',
          'hn',
          'hd',
          'hk',
          'pd',
          'pk',
          'z '];
% % %Adrian's model
% % %k_t d_t
% AA=[0 0;            %(3)
%     0 PdBar*DBar;        %(4)
%     PkBar*KBar 0;   %(5)
%     0 0;            %(6)
%     0 0;            %(7)
%     0 0;            %(8)
%     0 0;            %(9)
%     0 0;            %(10)
%     0 0;            %(11)
%     0 0;            %(12)
%     0 0;            %(13)
%     -KBar 0;        %(14)
%     0 -1;        %(15)
%     0 0            %(16)
%     0 1;   %(18)
%     1 0;  %(19)
%     ];
% %k_t-1, d_t-1
% BB=[0 0;            %(3)
%     0 0;        %(4)
%     0 0;   %(5)
%     0 0;            %(6)
%     0 0;            %(7)
%     0 0;            %(8)
%     0 0;            %(9)
%     0 0;            %(10)
%     0 0;            %(11)
%     0 0;            %(12)
%     -KBar 0;            %(13)
%     KBar*(1-delta) 0;        %(14)
%     0 (1-delta);        %(15)
%     0 0;            %(16) 
%     0 -(1-delta); %(18)
%     -(1-delta) 0]; %(19)
%     
% %c_n, c_d, i, y, r, w, k_n, k_d, k_k, h_n, h_d, h_k, p_d, p_k
% CC= [CnBar CdBar IBar -YBar 0 0 0 0 0 0 0 0 0 0; %(3)
%  0  -alpha*CdBar 0 0 KddBar*RBar 0 0 KddBar*RBar 0 0 0 0 0 0; %(4)
%  0 0 -alpha*IBar 0 KkdBar*RBar 0 0 0 KkdBar*RBar 0 0 0 0 0; %(5)
%  alpha*CnBar*RBar 0 0 0 alpha*CnBar*RBar 0 -KndBar 0 0 0 0 0 0 0; %(6)
%  -(1-alpha)*CnBar 0 0 0 0 WBar*HndBar 0 0 0 WBar*HndBar 0 0 0 0; %(7)
%  0 -alpha*CdBar 0 0 RBar*KddBar 0 0 RBar*KddBar 0 0 0 0 0 0; %(8)
%  0 -(1-alpha)*CdBar 0 0 0 WBar*HddBar 0 0 0 0 WBar*HndBar 0 0 0; %(9)
%  0 0 -alpha*IBar 0 RBar*KkdBar 0 0 0 RBar*KkdBar 0 0 0 0 0; %(10)
%  0 0 (1-alpha)*IBar 0 0 WBar*HkdBar 0 0 0 0 0 WBar*HkdBar 0 0; %(11)
%  0 0 0 0 0 0 0 0 0 HndBar HddBar HkdBar 0 0; %(12)
%  0 0 0 0 0 0 KndBar KddBar KkdBar 0 0 0 0 0; %(13)
%  0 0 0 0 0 0 0 0 alpha*KkdBar^alpha*HkdBar^(1-alpha) 0 0 (1-alpha)*KkdBar^alpha*HkdBar^(1-alpha) 0 0; %(14)
%  0 0 0 0 0 0 0 alpha*KddBar^alpha*HddBar^(1-alpha) 0 0 (1-alpha)*KddBar^alpha*HddBar^(1-alpha) 0 0 0; %(15)
%  -1 0 0 0 0 0 alpha 0 0 1-alpha 0 0 0 0; %(16) 
%  0 -delta 0 0 0 0 0 0 0 0 0 0 delta 0;
%  0 0 -delta 0 0 0 0 0 0 0 0 0 0 delta
%     ];
% %z_t
% DD=[0; %(1)
%     0; %(1)
%     0; %(1)
%     1; %(1)
%     1;
%     1;
%     1;
%     1;
%     1;
%     0;
%     0;
%     0;
%     0;
%     1;
%     ]

%Fiona's model
% Declaring the matrices
%k_t d_t
AA=[0 0;            %(3)
    0 PdBar*DBar;        %(4)
    PkBar*KBar 0;   %(5)
    0 0;            %(6)
    0 0;            %(7)
    0 0;            %(8)
    0 0;            %(9)
    0 0;            %(10)
    0 0;            %(11)
    0 0;            %(12)
    0 0;            %(13)
    -1 0;        %(14)
    0 -1;        %(15)
    0 0;            %(16)
    0 1;   %(18)
    1 0;  %(19)
    PkBar*KBar PdBar*DBar  %(20)
    ];
%k_t-1, d_t-1
BB=[0 0;            %(3)
    0 -PdBar*DBar*(1-delta);        %(4)
    -PkBar*KBar*(1-delta) 0;   %(5)
    0 0;            %(6)
    0 0;            %(7)
    0 0;            %(8)
    0 0;            %(9)
    0 0;            %(10)
    0 0;            %(11)
    0 0;            %(12)
    -KBar 0;            %(13)
    (1-delta) 0;        %(14)
    0 (1-delta);        %(15)
    0 0;            %(16) 
    0 -(1-delta); %(18)
    -(1-delta) 0 %(19)
    -(1-delta)*PkBar*KBar-RBar*KBar -(1-delta)*PdBar*DBar  %(20)
    ];
%c_n, c_d, i, y, r, w, k_n, k_d, k_k, h_n, h_d, h_k, p_d, p_k
CC= [CnBar CdBar IBar -YBar 0 0 0 0 0 0 0 0 0 0; %(3)
 0  -CdBar   0 0 0 0 0 0 0 0 0 0  PdBar*DBar*delta 0; %(4)
 0 0 -IBar 0 0 0 0 0 0 0 0 0 0  PkBar*KBar*delta; %(5)
 0 0 0 0 -1 0 alpha-1 0 0 1-alpha 0 0 0 0; %(6)
 0 0 0 0 0 -1 alpha 0 0 -alpha 0 0 0 0; %(7)
 0 0 0 0 -1 0 0 alpha-1 0 0 1-alpha 0 1 0; %(8)
 0 0 0 0 0 -1 0 alpha 0 0 -alpha 0 1 0; %(9)
 0 0 0 0 -1 0 0 0 alpha-1 0 0 1-alpha 0 1; %(10)
 0 0 0 0 0 -1 0 0 alpha 0 0 -alpha 0 1; %(11)
 0 0 0 0 0 0 0 0 0 HndBar HddBar HkdBar 0 0; %(12)
 0 0 0 0 0 0 KndBar KddBar KkdBar 0 0 0 0 0; %(13)
 0 0 0 0 0 0 0 0 alpha 0 0 1-alpha 0 0; %(14)
 0 0 0 0 0 0 0 alpha 0 0 1-alpha 0 0 0; %(15)
 -1 0 0 0 0 0 alpha 0 0 1-alpha 0 0 0 0; %(16) 
 0 -delta 0 0 0 0 0 0 0 0 0 0 delta 0; %(18)
 0 0 -delta 0 0 0 0 0 0 0 0 0 0 delta %(19)
 CnBar 0 0 0 -RBar*KBar -WBar 0 0 0 0 0 0 PdBar*DBar*delta PkBar*KBar*delta
    ];
%z_t
DD=[0; %(3)
    0; %(4)
    0; %(5)
    1; %(6)
    1; %(7)
    1; %(8)
    1; %(9)
    1; %(10)
    1; %(11)
    0; %(12)
    0; %(13)
    1; %(14)
    1; %(15)
    1; %(16)
    0; %(18)
    0; %(19)
    0; %(20)
    ];
%k_t+1, d_t+1
FF=[0 0; %(1)
    0 0]; %(2)
%k_t-1, d_t-1
HH=[0 0; %(1)
    0 0]; %(2)
%k_t, d_t
GG=[0 -beta*A/DBar; %(1)
    0 0]; %(2)
%y_t+1
JJ=[-beta*PdBar*(1-delta)/CnBar 0 0 0 0 0 0 0 0 0 0 0 beta*PdBar*(1-delta)/CnBar 0; %(1)
    -beta*PkBar*(1-delta)-beta*RBar 0 0 0 beta*RBar 0 0 0 0 0 0 0 0  beta*PkBar*(1-delta)]; %(2)
%y_t
KK=[PdBar/CnBar 0 0 0 0 0 0 0 0 0 0 0 -PdBar/CnBar 0; %(1)
    PkBar 0 0 0 0 0 0 0 0 0 0 0 0 -PkBar]; %(2)

LL=[0;0];
MM=[0;0];
NN=[rho];
Sigma=[sigma_e^2];
% Setting the options: 
[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);
HORIZON     = 32; % how far out should the impulse responses be calculated
PERIOD      = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
DO_PLOTS    = 0;  % if impulse response plots should be made, = 0, if not.
IMP_SELECT  = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
HP_SELECT   = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
GNP_INDEX   = 6; % Index of output among the variables selected for HP filter
% Calculation
do_it;
%% Plot

x = zeros(2,2000);
% yy(t) = [c(t); y(t); h(t); i(t); r(t); w(t)]
yy = zeros(14,2000);
% zz(t) = [z(t)]
zz = zeros(1,2000);

for t = 2:2000
    zz(t) = rho*zz(t-1) + randn*sigma_e;
    x(:,t) = PP * x(:,t-1) + QQ * zz(t);
    yy(:,t) = RR * x(:,t-1) + SS * zz(t);
end

% plot cyclical component of y vs. time
figure
plot(501:2000,yy(4,501:2000));
title('cyclical component output');
xlabel('time');
ylabel('output.');

wanted=yy(1:4,501:2000)';

%The order is different from the data one, data is y, cn, cd, i, and the
%model table's order is cn cd i y.
CorrMatrix_model = corr(wanted);
RelVol_model = std(wanted)/std(wanted(:,4));
table_model=[RelVol; CorrMatrix];

%The model doesn't capture the correlation very well. It over estimate the
%economy's consumption and underestimate the volitility of investments.
%It captures 65% of the relative volatility bewteeen non-durable goods and
%durable goods. For specifics please see the table. 