% Solve a set of ordinary differential equations using the most simple
% numerical method, i.e. the Euler method. The method is described at:
% http://en.wikipedia.org/wiki/Euler_method
%
% Function's input and output paramaeters should conform to the general integ
% "interface" as described in README.
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


function output = integ_euler(model, initial_condition, t_start, t_stop, t_step, outputf, param)


% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);
output = [t_start; initval];

s = initial_condition;
for t = t_start : t_step : t_stop-t_step,
    sp = feval(model, s, t, param);
    s = s + sp * t_step;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
    
end % for

end % function
