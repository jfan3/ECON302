% SOL_OUT.M prints the coefficients of the decision rules,
% delivered by SOLVE.M.
% It is assumed, that VARNAMES, a matrix with m+n+k rows has
% been set, containing the names of all the variables.
% This program overwrites m_states, k_exog and n_endog.

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
[m_states,k_exog] = size(QQ);
[n_endog,k_exog] = size(SS);
disp('Exogenous states z(t):');
disp(VARNAMES((m_states+n_endog+1):(m_states+n_endog+k_exog),:));
disp(' ');
disp('Endogenous states x(t):');
disp(VARNAMES(1:m_states,:));
disp(' ');
disp('PP: Recursive equilibrium law of motion for x(t) on x(t-1):');
disp(PP);
disp('QQ: Recursive equilibrium law of motion for x(t) on z(t):');
disp(QQ);
disp(' ');
disp('Other endogenous variables y(t):');
disp(VARNAMES((m_states+1):(m_states+n_endog),:));
disp(' ');
disp('RR: Recursive equilibrium law of motion for y(t) on x(t-1):');
disp(RR);
disp('SS: Recursive equilibrium law of motion for y(t) on z(t):');
disp(SS);
disp(' ');
