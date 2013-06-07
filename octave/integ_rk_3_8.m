% Solve a set of ordinary differential equations using the Runge - Kutta 3/8 method.
% The method is described at:
% http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
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

function output = integ_rk_3_8(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Elements of the Butcher tableau:
A = [0, 0, 0, 0; 1, 0, 0, 0; -1, 3, 0, 0; 3, -3, 3, 0] / 3;
B = [0.125, 0.375, 0.375, 0.125];
C = [0, 1, 2, 3]' / 3;

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C);

end % function  
