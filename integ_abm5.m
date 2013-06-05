% Solve a set of ordinary differential equations using the
% 5-step Adams - Bashforth - Moulton predictor - corrector method. The method is described at:
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
% and
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods
% Coefficients for the 6-step Adams - Moulton method are available at:
% http://www.cs.rpi.edu/~flaherje/pdf/ode5.pdf
%
% As the method is not self starting, the first points are calculated by the 
% 4th order Runge - Kutta method.
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


function output = integ_abm5(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Coefficients for the 5-step Adams - Bashforth method:
ab_coef = [1901, -2774, 2616, -1274, 251] / 720;

% Coefficients for the 6-step Adams - Moulton method:
am_coef = [475, 1427, -798, 482, -173, 27] / 1440;

output = aux_abm_general(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, ab_coef, am_coef);

end % function  
