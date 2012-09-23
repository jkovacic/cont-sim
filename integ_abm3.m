% Solve a set of ordinary differential equations using the
% 3-step Adams - Bashforth - Moulton predictor - corrector method. The method is described at:
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


function output = integ_abm3(model, initial_condition, t_start, t_stop, t_step, outputf, param)


% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% This method is not self starting at the first N points must be calculated
% by another method. The 4th order Runge - Kutta method is chosen.
% See integ_rk4.m for more details about the method.

% This design allows easier reusability in other non selfstarting methods.

% number of initial points to be calculated using a self starting method:
N = 3;

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);
output = [t_start; initval];

% The 4th order Runge - Kutta method to solve the first N points
s = initial_condition;
upper = t_start + (N-1)*t_step;

% Make sure that not too many points will be calculated if t_step is too large
if (upper > (t_stop-t_step) )
    upper = t_stop-t_step;
end %if

sd1 = 0;

for  t = t_start : t_step : upper
    % used by the Adams - Bashforth method as "previous" points:
    sd2 = sd1;
    sd1 = feval(model, s, t, param);
    
    k1 = sd1 * t_step;
    k2 = feval(model, s+0.5*k1, t+0.5*t_step, param) * t_step;
    k3 = feval(model, s+0.5*k2, t+0.5*t_step, param) * t_step;
    k4 = feval(model, s+k3, t+t_step, param) * t_step;

    s = s + (k1 + 2*k2 + 2*k3 + k4) / 6;

    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
end %for


% Now the Adams - Bashforth - Moulton method can start

for t = upper+t_step : t_step : t_stop-t_step,

    % 3-step Adams - Bashforth method is used as a predictor:
    sd = feval(model, s, t, param);
    p = s + t_step * (23*sd - 16*sd1 + 5*sd2) / 12;
    
    % and corrected by the 3-step Adams - Moulton method:
    pd = feval(model, p, t+t_step, param);
    s = s + t_step * (9*pd + 19*sd - 5*sd1 + sd2) / 24;

    sd2 = sd1;
    sd1 = sd;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
    
end % for

end % function 
 
