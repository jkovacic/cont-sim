% Solve a set of ordinary differential equations using the most simple
% numerical method, i.e. the Euler method. The method is described at:
% http://en.wikipedia.org/wiki/Euler_method
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


function output = integ_euler(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)


% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% The first set of output values at t = t_start:
ut = feval(inputf, t_start);
initval = feval(outputf, initial_condition, ut, t_start, param);

% To improve efficiency, preallocate the buffer for output:
[OROWS, IGNORED] = size(initval);
output = aux_preallocate_output(OROWS+1, t_start, t_stop, t_step);

% Fill the 'output' with the initial value:
output(:, 1) = [t_start; initval];

s = initial_condition;
% Current index within 'output'
idx = 2;
for t = t_start : t_step : t_stop-t_step,
    sd = feval(model, s, ut, t, param);
    s = s + sd * t_step;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calculating output values and applying their time stamp.
    ut = feval(inputf, t+t_step);
    val = feval(outputf, s, ut, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    
    % update 'idx'
    idx = idx+1;
end % for

end % function
