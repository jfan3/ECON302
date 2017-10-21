% README.M contains hints on how to use the programs in this directory.
% Start MATLAB and type
% readme
% on a single line.  Alternatively, inspect this file with any
% text editor.

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.

disp('The files in this directory perform the calculations described');
disp('in Harald Uhlig (1995), "A Toolkit for Analyzing Nonlinear Dynamic');
disp('Stochastic Models Easily", Discussion Paper, Institute for');
disp('Empirical Macroeconomics, Federal Reserve Bank of Minneapolis.');
disp('To see quickly, how they work, start MATLAB and type');
disp('example1');
disp('to calculate through example 1, i.e. Hansen RBC model, or type');
disp('example2');
disp('to calculate through example 2, i.e. the Lettau-Uhlig RBC model');
disp('with habit formation.  Use the files example1.m or example2.m');
disp('as templates for your own work.  Alternatively, declare all');
disp('needed matrices and type in');
disp('do_it');
disp('to do all calculations.');
disp(' ');
disp('The files in this directory are:');
disp('do_it.m      : does it all, once all needed matrices are defined.');
disp('example1.m   : Hansen RBC model, example 1 in the paper.');
disp('example2.m   : habit-RBC model, example 2 in the paper.');
disp('impresp.m    : calculates and shows impulse responses to shocks.');
disp('mom_out.m    : produces output.  To be called after moments.m');
disp('moments.m    : calculates second moment properties.');
disp('readme.m     : this file.  It tells you what to do.');
disp('sol_out.m    : produces output.  To be called after solve.m');
disp('solve.m      : solves for the recursive equilibrium law of motion.');
disp(' ');
disp('All files are extensively documented.  Type');
disp('help filename');
disp('in MATLAB to get more information. Note that these files');
disp('set some additional variables, which you may have used before:');
disp('thus, be careful not to use names appearing in the programs.');
disp('If you discover a serious mistake, please contact me at uhlig@kub.nl.');
disp('If you have a question, please do NOT contact me.  Please read the paper.');
disp('Feel free to copy, modify and use these files at your own risk.');
disp('There is absolutely no guarantee that this stuff works the way it is');
disp('supposed to.  Have fun.');
disp('                Harald Uhlig, Minneapolis, June 1995');
