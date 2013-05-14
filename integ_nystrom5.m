% Solve a set of ordinary differential equations using the 5th order Nystroem's method
% The method is described at:
% http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_nystrom.c
%
% Function's input and output paramaeters should conform to the general integ
% "interface" as described at https://github.com/jkovacic/cont-sim/wiki/Basic-instructions.
%
% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps

function output = integ_nystrom5(model, initial_condition, t_start, t_stop, t_step, outputf, param)



% Elements of the Butcher tableau:
A = [ 0, 0, 0, 0, 0, 0; ...
      1/3, 0, 0, 0, 0, 0; ...
      0.16, 0.24, 0, 0, 0, 0; ...
      0.25, -3, 15/4, 0, 0, 0; ...
      2/27, 10/9, -50/81, 8/81, 0, 0; ...
      0.08, 0.48, 2/15, 8/75, 0, 0];
B = [23, 0, 125, 0, -81, 125] / 192;
C = [0, 1/3, 0.4, 1, 2/3, 0.8]';

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, outputf, param, A, B, C);

end % function  
