% MOM_OUT produces output from the calculations done with MOMENTS.M,
% which is assumed to have been run just before.
% This program should be modified to suit tastes and needs.  Some options
% are given in the first few lines of this program.
% It is assumed that
% VARNAMES, a matrix with (m+n+k) rows, containing the variable names, has been set.
% The program overwrites freq1, diag_select, hndl, var_index

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.

% Options:
DO_HP_GRAPH = 0; % Set to = 1 to see plots of the spectral densities.
DO_DISP1    = 1; % Set to = 1 to see printout of the autocorrelation matrix. 
DO_DISP2    = 0; % Set to = 1 to see printout of the variance-covariance matrix.
DO_DISP3    = 1; % Set to = 1 to see printout of the vector of variances.

if DO_HP_GRAPH,
   freq1 = N_GRIDPOINTS/24;
   diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
   hndl = plot(freqs(freq1:N_GRIDPOINTS/2), real(svv_raw(freq1:N_GRIDPOINTS/2,diag_select)));
   set(hndl,'LineWidth',2);
   for var_index = 1:n_select,
      text(freqs(freq1), real(svv_raw(freq1,diag_select(var_index))),...
           VARNAMES(HP_SELECT(var_index),:));
   end;
   title('Spectral densities, unfiltered');
   xlabel('Frequency');
   disp('Inspect figure. Hit key when ready...');
   pause;
diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
   hndl = plot(freqs(1:N_GRIDPOINTS/2), real(svv_fil(1:N_GRIDPOINTS/2,diag_select)));
   set(hndl,'LineWidth',2);
   for var_index = 1:n_select,
      text(freqs(N_GRIDPOINTS/20), real(svv_fil(N_GRIDPOINTS/20,diag_select(var_index))),...
           VARNAMES(HP_SELECT(var_index),:));
   end;
   title('Spectral densities, Hodrick-Prescott filtered');
   xlabel('Frequency');
   disp('Inspect figure. Hit key when ready...');
   pause;
end;
disp(' ');
disp('The variables are:');
disp(VARNAMES(HP_SELECT,:));
disp(' ');
if DO_DISP1,
   disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
   for var_index = 1 : n_select,
      disp(sprintf('  %5.2f',autcor_fil(var_index,:)));
   end; 
   disp(sprintf('  %5.0f',autcor_fil(n_select+1,:)));
   disp(' ');
end;
if DO_DISP2,
   disp('Variance-Covariance Matrix, HP-filtered series:');
   for var_index = 1 : n_select,
      disp(sprintf(' %6.3f',covmat_fil(var_index,:)));
   end;
   disp(' ');
end;
if DO_DISP3,
   disp('Standard deviations, HP-filtered series:');
   disp(sqrt(varvec_fil));
   disp(' ');
   disp('Relative volatilities, HP-filtered series:');
   disp(sqrt(varvec_fil)/sqrt(varvec_fil(GNP_INDEX)));
end;
  
