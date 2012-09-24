% Solve a set of ordinary differential equations using the
% 2-step Adams - Bashforth - Moulton predictor - corrector method. The method is described at:
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
% and
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods
%
% As the method is not self starting, the first points are calculated by the 
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
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps


function output = integ_abm2(model, initial_condition, t_start, t_stop, t_step, outputf, param)

% Coefficients for the 2-step Adams - Bashforth method:
ab_coef = [1.5, -0.5];

% Coefficients for the 3-step Adams - Moulton method:
am_coef = [5, 8, -1] / 12;

output = aux_abm_general(model, initial_condition, t_start, t_stop, t_step, outputf, param, ab_coef, am_coef);

end % function 
