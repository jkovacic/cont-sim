% Solve a set of ordinary differential equations using the most common
% numerical method, i.e. the 4th order Runge - Kutta method. The method is described at:
% http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Common_fourth-order_Runge.E2.80.93Kutta_method
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

function output = integ_rk4(model, initial_condition, t_start, t_stop, t_step, outputf, param)


% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);
output = [t_start; initval];

s = initial_condition;
for t = t_start : t_step : t_stop-t_step
    k1 = feval(model, s, t, param) * t_step;
    k2 = feval(model, s+0.5*k1, t+0.5*t_step, param) * t_step;
    k3 = feval(model, s+0.5*k2, t+0.5*t_step, param) * t_step;
    k4 = feval(model, s+k3, t+t_step, param) * t_step;
    
    s = s + (k1 + 2*k2 + 2*k3 + k4) / 6;

    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
end % for

end % function 
