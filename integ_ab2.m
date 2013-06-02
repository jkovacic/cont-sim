% Solve a set of ordinary differential equations using the
% 2-step Adams - Bashforth method. The method is described at:
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
%
% As the method is not sel starting, the first points are calculated by the 
% 4th order Runge - Kutta method.
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
%   inputf - name of the function that returns the external input at the specified time
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps


function output = integ_ab2(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Coefficients for the 2-step Adams - Bashforth method:
coef = [1.5, -0.5];

output = aux_ab_general(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, coef);

end % function
