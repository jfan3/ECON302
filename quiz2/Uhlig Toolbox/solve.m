% SOLVE.M solves for the decision rules in a linear system,
% which is assumed to be of the form
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,
% where it is assumed that x(t) is the endogenous state vector,
% y(t) the other endogenous variables and z(t) the exogenous state
% vector.  It is assumed that the row dimension of AA is at least as large as
% the dimensionality of the endogenous state vector x(t).  
% The program solves for the equilibrium law of motion
% x(t) = PP x(t-1) + QQ z(t)
% y(t) = RR x(t-1) + SS z(t).
% To use this program, define the matrices AA, BB, .., NN.
% SOLVE.M then calculates PP, QQ, RR and SS.  It also calculates
% WW with the property [x(t)',y(t)',z(t)']=WW [x(t)',z(t)'].
% 
% A few additional variables are used
% (TOL, l_equ, m_states, n_endog, k_exog,
% CC_plus, CC_0, Psi_mat, Gamma_mat, Theta_mat, Xi_mat,
% Xi_eigvec, Xi_eigval, Xi_eigsort, Xi_sortindex,
% Xi_sortsub, entry1, last_in, last2_in, entry2,
% Omega_piece, Lambda_mat, Omega_mat, VV,
% LLNN_plus_MM, QQSS_vec), overwriting variables with the same names
% that might have been used before.
% Source: H. Uhlig (1995) "A Toolkit for Solving Nonlinear Dynamic
% Stochastic Models Easily," Discussion Paper, Institute for
% Empirical Macroeconomis, Federal Reserve Bank of Minneapolis

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.

% OPTIONS:
TOL = .000001; % Roots smaller than TOL are regarded as zero.

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
k_exog = min(size(NN));
if rank(CC)<n_endog,
 disp('SOLVE.M: Sorry!  Rank(CC) needs to be at least n! Cannot solve for PP.');
else
 CC_plus = pinv(CC);
 CC_0 = (null(CC'))';
 Psi_mat = [ CC_0*AA
            FF-JJ*CC_plus*AA];
 if rank(Psi_mat)<m_states,
  disp('SOLVE.M: Sorry!  Psi singular! Cannot solve for PP.');
 else
  Gamma_mat = Psi_mat \ [ - CC_0*BB
                          JJ*CC_plus*BB-GG+KK*CC_plus*AA ];
  Theta_mat = Psi_mat \ [ zeros(l_equ-n_endog,m_states)
                          KK*CC_plus*BB - HH                ];
  Xi_mat = [ Gamma_mat,     Theta_mat
             eye(m_states), zeros(m_states) ];
  [Xi_eigvec,Xi_eigval] = eig(Xi_mat);
  if rank(Xi_eigvec)<m_states,
   disp('SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.');
  else
   [Xi_eigsort,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
   entry1 = 1; % index in Xi_sortindex of the first root to be used
   % The next loop eliminates spurious zero roots, which
   % are created by postmultiplying the equation C0*A*P+C0*B = 0 with P,
   % which needed to be done to obtain a matrix quadratic equation which
   % is solvable as a standard eigenvalue problem.
   % The postmultiplication is problematic only if P is singular, i.e.
   % if P has zero roots. Note that there should not be zero roots anyways,
   % if the state space has been chosen to be of minimal dimension.
   % Thus, the following loop is not problematic.
   while ( (entry1<=l_equ-n_endog)&(Xi_eigsort(entry1)<TOL))
    entry1 = entry1 + 1;
   end;
   if ((m_states > 1) & (entry1 <= m_states)),
    last_in = Xi_sortindex(entry1 + m_states-1);  % rowindex in Xi_eigval of largest root to be used
    last2_in = Xi_sortindex(entry1 + m_states-2);
    if imag(Xi_eigval(last_in,last_in))~=0,
     if Xi_eigval(last2_in,last2_in)~=...
           conj(Xi_eigval(last_in,last_in)),
      disp('SOLVE.M: I will drop the lowest eigenvalue to get real PP. Hope that is ok.');
      entry = entry1+1;
     end;
    end;
   end;
   if Xi_eigsort(entry1+m_states-1) > 1,
    disp('SOLVE.M: You may be in trouble.  There are not enough stable roots for PP.');
   end;
   if entry1 <= m_states,
    if Xi_eigsort(entry1+m_states) < 1,
     disp('SOLVE.M: Too many stable roots. I will use the smallest roots. Hope that is ok.');
    end;
   end;
   entry2=0;
   Omega_piece = [];
   if l_equ > n_endog,
    while ((Xi_eigsort(entry1+entry2)<TOL)&(entry2<m_states)),
       entry2 = entry2+1;
    end;
    % The following piece can only arise if the endogenous state
    % space was not chosen to be of minimal size.
    if entry2 > 0,
     if entry2 >  m_states - 1,
      disp('SOLVE.M: WARNING! WARNING! WARNING: TOO MANY ZERO ROOTS.  THUS,');
      disp('SOLVE.M:     I PROBABLY CANNOT FIND A SENSIBLE SOLUTION.  PLEASE');
      disp('SOLVE.M:     ELIMINATE UNNECCESSARY ENDOGENOUS STATES.  I WILL');
      disp('SOLVE.M:     PROCEED BUT P PROBABLY DOES NOT SATISFY ITS EQUATION!');
      entry2 = 0;
     else
      % the linear part of the matrix quadratic equation for P gives the
      % restriction that the null space of PP must be a subset of the
      % nullspace of CC_0*BB.  Thus, to find the null space of PP,
      % intersect the nullspace of CC_0*BB with the
      % null space of Xi_mat, projected on its lower half:
      Omega_piece = null([[0*CC_0*BB,CC_0*BB];Xi_mat]);
      Omega_piece = [ zeros(m_states), eye(m_states) ]*Omega_piece; % chop to lower part
      Omega_piece = orth(Omega_piece); % To make sure that there are no linear dependencies
      % This works fine, except if eigenvalues identified to be zero via
      % the loop for entry2 above are actually not regarded as zero by
      % null(.). In that case, I guess that the tolerance TOL was set
      % too high and thus entry2 is too high. Reset entry2 accordingly:
      entry2 = min(size(Omega_piece));
     end;
    end;  
   end;  
   Xi_sortsub = Xi_sortindex((entry1+entry2):(entry1+m_states-1));
   Lambda_mat = Xi_eigval(Xi_sortsub,Xi_sortsub);
   Lambda_mat = [zeros(entry2,entry2),          zeros(entry2,m_states-entry2)
                 zeros(m_states-entry2,entry2), Lambda_mat                    ];
   Omega_mat  = [Omega_piece, Xi_eigvec((m_states+1):(2*m_states),Xi_sortsub)];
   if rank(Omega_mat)<m_states,
    disp('SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.');
   else
    PP = Omega_mat*Lambda_mat/Omega_mat;
    PP_imag = imag(PP);
    PP = real(PP);
    if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
     disp('SOLVE.M: PP is complex.  I proceed with the real part only.  Hope that is ok.');
    end;
    RR = - CC_plus*(AA*PP+BB);
    VV = [ kron(eye(k_exog),AA),   kron(eye(k_exog),CC)
           kron(NN',FF)+kron(eye(k_exog),(FF*PP+JJ*RR+GG)), ...
                                   kron(NN',JJ)+kron(eye(k_exog),KK) ];
    if rank(VV) < k_exog*(m_states+n_endog),
     disp('SOLVE.M: Sorry! V is not invertible.  Cannot solve for QQ and SS.');
    else
     LLNN_plus_MM = LL*NN + MM;
     QQSS_vec = - VV \ [ DD(:)
                         LLNN_plus_MM(:) ];
     QQ = reshape(QQSS_vec(1:m_states*k_exog),m_states,k_exog);
     SS = reshape(QQSS_vec((m_states*k_exog+1):((m_states+n_endog)*k_exog)),n_endog,k_exog);
     WW = [ eye(m_states)         , zeros(m_states,k_exog)
            RR*pinv(PP)           , (SS-RR*pinv(PP)*QQ) 
            zeros(k_exog,m_states), eye(k_exog)            ];
    end;
   end;
  end;
 end;
end;


       
     
      

    
  
      
      

     
