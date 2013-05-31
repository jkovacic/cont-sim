% Solve a set of ordinary differential equations using the
% 1-step Adams - Bashforth - Moulton predictor - corrector method (actually the 
% corrected Euler method). The method is described at:
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
% and
% http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods
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


function output = integ_abm1(model, initial_condition, t_start, t_stop, t_step, outputf, param)


% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);

% To improve efficiency, preallocate the buffer for output:
[OROWS, IGNORED] = size(initval);
output = aux_preallocate_output(OROWS+1, t_start, t_stop, t_step);
% Fill the 'output' with initial values
output(:, 1) = [t_start; initval];

s = initial_condition;

% Current index within 'output'
idx = 2;
for t = t_start : t_step : t_stop-t_step,

    % 1-step Adams - Bashforth method (actually the Euler method) is used as a predictor:
    sd = feval(model, s, t, param);
    p = s + t_step * sd;
    
    % and corrected by the 1-step Adams - Moulton method:
    pd = feval(model, p, t+t_step, param);
    s = s + t_step * pd;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    
    % update 'idx'
    idx = idx+1;
end % for

end % function 
 
 
