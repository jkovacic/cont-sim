% Solve a set of ordinary differential equations using the 2nd order Ralston's method,
% also known as the  alternative 2nd order Heun method 
% (similar to the 2nd Runge - Kutta method, alpha=2/3)
% The method is described at:
% http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_2.html
% and http://en.wikipedia.org/wiki/Heun%27s_method
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

function output = integ_ralston2(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Elements of the Butcher tableau:
[A, B, C] = aux_rk2_butcher(2/3);

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C);

end % function   
