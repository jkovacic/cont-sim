% Solve a set of ordinary differential equations using the 6th order Butcher's method.
% The method is described at:
% http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_butcher.c
%
% Function's input and output parameters should conform to the general integ
% "interface" as described at https://github.com/jkovacic/cont-sim/wiki/Basic-instructions.
%
% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step
%   inputf - name of the function that returns the external input at the specified time
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps

function output = integ_butcher6(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Elements of the Butcher tableau:
A = [ 0, 0, 0, 0, 0, 0, 0; ... 
      1/3, 0, 0, 0, 0, 0, 0; ... 
      0, 2/3, 0, 0, 0, 0, 0; ...
      1/12, 1/3, -1/12, 0, 0, 0, 0; ...
      -0.0625,  1.125, -0.1875, -0.375, 0, 0, 0; ...
      0, 1.125, -0.375, -0.75, 0.5, 0, 0; ...
      9/44, -9/11, 63/44, 18/11, -16/11, 0, 0 ];
B = [11, 0, 81, 81, -32, -32, 11] / 120;
C = [0, 1/3, 2/3, 1/3, 0.5, 0.5, 1]';

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C);

end % function  
