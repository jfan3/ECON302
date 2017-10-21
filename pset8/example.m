% EXAMPLE. calculates the habit-formation model of Lettau and Uhlig (1995), 
% "Can habit formation be reconciled with business cycle facts?"

disp('EXAMPLE2: A real business cycle model with Campbell-Cochrane-type');
disp('          habit formation, see Lettau, M. and Uhlig, H. , ');
disp('          "Can Habit Formation be Reconciled With Business Cycle Facts?",');
disp('          draft, CentER, Tilburg University (1995).');

% Setting parameters:

N_bar     = 1.0/3; % Steady state hours worked is one third of total time endowment
Z_bar     = 1;     % Normalization
rho       = .36;   % Capital share
delta     = .025;  % Depreciation rate
R_bar     = 1.01;  % One percent real interest per quarter
eta       = 2.372; % The value recommended by Campbell and Cochrane (1994),"By force of habit..."
phi       = 0.97;  % The value recommended by Campbell and Cochrane (1994),"By force of habit..."
psi       = .95;   % First order autocorrelation of the technology parameter z(t)
sigma_eps = 0.763; % Source: Prescott (1986), Theory ahead of measurement, page 16.  Units: percent
S_bar     = 0.05;  % Steady state consumption surplus, recommended by Campbell and Cochrane (1994).

% Calculating steady state relationships:

betta   = 1.0/R_bar;
YK_bar  = (R_bar + delta - 1)/rho;  % = Y_bar / K_bar
K_bar   = (YK_bar / Z_bar)^(1.0/(rho-1)) * N_bar;
Y_bar   = YK_bar * K_bar;
C_bar   = Y_bar - delta*K_bar;
D_bar   = rho * YK_bar;
W_bar   = (1-rho) * Y_bar/N_bar;
A       = (S_bar * C_bar)^(-eta) * (1 - rho) * Y_bar/N_bar;
lambda  = 1.0/S_bar - 1.0; % Using a steady state relationship taken
   % from functional form in Campbell and Cochrane (1994),
   % see also Lettau and Uhlig (1995). In Uhlig (1995), 
   % "A toolkit...", no restriction on lambda is stated.
X_bar   = C_bar*(1 - S_bar); % Steady state habit


% Declaring the matrices. 

% Endogenous state variables "x(t)" is capital and q(t) = phi s(t) - lambda c(t)

VARNAMES = ['capital       ',
            'q(t)          ',
            'cons. surplus ',
            'consumption   ',
            'output        ',
            'labor         ',
            'dividends     ',
            'wages         ',
            'interest      ',
            'habit         ',
            'technology    '];

% Translating into coefficient matrices.  
% The equations are, conveniently ordered:
% 1) 0 = - K k(t) + (D+1-delta)K k(t-1) - C c(t) + DK d(t) + WN (w(t) + n(t))
% 2) 0 = rho k(t-1) - y(t) + (1-rho) n(t) + z(t)
% 3) 0  = q(t-1) - s(t) + lambda c(t)   --> Note the replacement with q(t)
% 4) 0 = - q(t) + phi s(t) - lambda c(t)  --> Extra equations to define q(t)
% 5) 0 = -k(t-1) + y(t) - d(t)
% 6) 0 = y(t) - n(t) - w(t)
% 7) 0 = - eta (s(t)+c(t)) + y(t) - n(t) 
% 8) 0 = - rho Y/K k(t-1) + rho Y/K y(t) - R r(t)
% 9) 0 = S_bar s(t) + X_bar/C_bar ( x(t) - c(t) )  
%            --> Loglinearizing definition of S(t) to get x(t)
% 10) 0 = E_t [ - eta (s(t+1) + c(t+1) )+ r(t+1) + eta (s(t) +  c(t)) ]
% 11) z(t+1) = psi z(t) + epsilon(t+1)
% CHECK: 11 Equations, 11 variables.


AA = [ - K_bar,  0   % Equ.1)
             0,  0   % Equ.2)
             0,  0   % Equ.3)
             0, -1   % Equ.4)
             0,  0   % Equ.5)
             0,  0   % Equ.6)
             0,  0   % Equ.7)
             0,  0   % Equ.8)
             0,  0 ];% Equ.9)

BB = [ (D_bar+1-delta)*K_bar,  0    % Equ.1)
                         rho,  0    % Equ.2)
                           0,  1    % Equ.3)
                           0,  0    % Equ.4)
                          -1,  0    % Equ.5)
                           0,  0    % Equ.6)
                           0,  0    % Equ.7)
                - rho*YK_bar,  0    % Equ.8)
                           0,  0 ]; % Equ.9)


% Ordering of other endogenous variables:
%      s,       c,      y,       n,       d,       w,      r,   x

CC = [ 0,  -C_bar,   0,W_bar*N_bar,D_bar*K_bar,W_bar*N_bar,0,   0   % Equ.1)
       0,       0,     -1, (1-rho),       0,       0,      0,   0   % Equ.2)
      -1,  lambda,      0,       0,       0,       0,      0,   0   % Equ.3)
     phi, -lambda,      0,       0,       0,       0,      0,   0   % Equ.4)
       0,       0,      1,       0,      -1,       0,      0,   0   % Equ.5)
       0,       0,      1,      -1,       0,      -1,      0,   0   % Equ.6)
    -eta,    -eta,      1,      -1,       0,       0,      0,   0   % Equ.7)
       0,       0,rho*YK_bar,    0,       0,       0, -R_bar,   0   % Equ.8)
   S_bar,-X_bar/C_bar,  0,       0,       0,       0,0,X_bar/C_bar];% Equ.9)

DD = [ 0
       1
       0
       0
       0
       0
       0
       0
       0 ];

FF = [ 0 0 ];

GG = [ 0 0 ];

HH = [ 0 0 ];

% Ordering of other endogenous variables:
%         s,    c, y,  n,  d,  w,  r,  x

JJ = [ -eta, -eta, 0,  0,  0,  0,  1,  0 ];

KK = [  eta,  eta, 0,  0,  0,  0,  0,  0 ];

LL = [ 0 ];

MM = [ 0 ];

NN = [ psi ];

Sigma = [ sigma_eps^2 ];

% Setting the options:
  
[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);
HORIZON     = 32; % how far out should the impulse responses be calculated
PERIOD      = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
DO_PLOTS    = 1;  % if impulse response plots should be made, = 0, if not.
IMP_SELECT  = 1:(m_states+n_endog+k_exog);
   %  a vector containing the indices of the variables to be plotted
HP_SELECT   = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
GNP_INDEX   = 5; % Index of output among the variables selected for HP filter


% Starting the calculations:

do_it;


 