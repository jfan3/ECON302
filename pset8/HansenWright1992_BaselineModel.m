% Uses Uhlig's MATLAB Toolkit to solve baseline RBC model from lecture.
% First, parameters are set and the steady state is calculated. Next, the
% matrices are declared. In the last line, the model is solved and
% analyzed by calling DO_IT.m

clear all clc

% Calibration and calculating the steady state:

h_bar     = 1/3;  % Steady state employment
r_bar     = 0.035; % steady state return on capital before depreciation
z_bar     = 1; % Normalization
theta     = 0.36; % Capital's share of total income
delta     = 0.025; % Depreciation rate for capital
rhoz      = 0.95; % autocorrelation of technology shock
sigma_eps = 0.7; % Standard deviation of technology shock.  Units: Percent.

bet     = (r_bar+1-delta)^(-1);  % Discount factor beta
ky_bar  = theta/r_bar; % steady state capital to output ratio
cy_bar  = 1-delta*ky_bar; % steady state consumption to output ratio
k_bar   = (r_bar/theta)^(1/(theta-1))*h_bar; % steady state capital
y_bar   = k_bar/ky_bar; % steady state output
c_bar   = cy_bar*y_bar; % steady state consumption
w_bar   = (1-theta)*k_bar^theta*h_bar^(-theta);
i_bar   = delta*k_bar;
A       = (1-h_bar)*(1-theta)/(h_bar*cy_bar); % preference parameter


% Declaring the matrices.

VARNAMES = ['capital    ',
    'consumption',
    'output     ',
    'labor      ',
    'investment ',
    'rental rate',
    'wage rate  ',
    'technology '];

% Endogenous state variables "x(t)": k(t)
% Endogenous other variables "y(t)": [c(t), y(t), h(t), i(t), r(t), w(t)]'
% Exogenous state variables  "z(t)": z(t)
% Find matrices for format:
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

AA = [   0
      -k_bar
         0
         0
         0
         0];

BB = [ r_bar*k_bar
       k_bar*(1-delta)
       theta
       0
       theta-1
       theta];

% Order:   c(t)        y(t)         h(t)           i(t)         r(t)         w(t)
CC = [  -c_bar,        0,       w_bar*h_bar,    -i_bar,     r_bar*k_bar,  w_bar*h_bar
          0,           0,           0,           i_bar,         0,            0
          0,          -1,       1-theta,          0,            0,            0
       A*c_bar,        0,       w_bar*h_bar,      0,            0,      -w_bar*(1-h_bar)
          0,           0,       (1-theta),        0,           -1,            0
          0,           0,         -theta,         0,            0,           -1];

DD = [ 0
       0
       1
       0
       1
       1];

FF = [ 0 ];

GG = [ 0 ];

HH = [ 0 ];

% Order:   c(t+1) y(t+1) h(t+1) i(t+1)   r(t+1)   w(t+1)
JJ =     [  -1,    0,     0,     0,   bet*r_bar,   0 ];

% Order:   c(t) y(t) h(t) i(t) r(t) w(t)
KK =     [ 1,   0,   0,   0,   0,   0 ];

LL = [ 0 ];

MM = [ 0 ];

NN = [rhoz];

Sigma = [ sigma_eps^2  ];

% Setting the options:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);
HORIZON    = 32; % how far out should the impulse responses be calculated
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
DO_PLOTS   = 1;  % if impulse response plots should be made, = 0, if not.
IMP_SELECT = 1:(m_states+n_endog+k_exog);
% a vector containing the indices of the variables to be plotted
IMP_SELECT = [1:8];  % variables in the impulse-response plots
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
GNP_INDEX  = 3; % Index of output among the variables selected for HP filter


% Starting the calculations:

do_it;



%% Plot

x = zeros(1,2000);
% yy(t) = [c(t); y(t); h(t); i(t); r(t); w(t)]
yy = zeros(6,2000);
% zz(t) = [z(t)]
zz = zeros(1,2000);

for t = 2:2000
    zz(t) = rhoz*zz(t-1) + randn*sigma_eps;
    x(t) = PP * x(t-1) + QQ * zz(t);
    yy(:,t) = RR * x(t-1) + SS * zz(t);
end

% plot h vs. w
figure
plot(yy(3,501:2000),yy(6,501:2000),'x');
title('Standard Model');
xlabel('Hours');
ylabel('Prod.');
axis([-10 10 -10 10]);