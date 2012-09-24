% A convenience function the first N solutions of a set of ordinary differential equations,
% typically used as a start of multi step methods (e.g. the Adams family).
% This function should not be called directly

% The function is similar to integ_rk4.m. However, this implementation is suited to
% multi step methods, for instance, it also maintains a buffer of the last N values
% od model(sm, tn), etc.

% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   upper - stop time of the simulation run
%   t_step - fixed time step
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%   N - number of previous states to be stored in H
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps
%   s - vector of states at the last point (at t = upper)
%   H - buffer of values of model(sn, tn) for the last N points

function [output, s, H] = aux_rk4(model, initial_condition, t_start, upper, t_step, outputf, param, N)

[STATE_ROWS, STATE_COLS] = size(initial_condition);

% A buffer to store previous N-1 values of model(sn, tn), 
% required by the Adams - Bashforth method
H = zeros(STATE_ROWS, (N-1) * STATE_COLS);

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);
output = [t_start; initval];

s = initial_condition;


for  t = t_start : t_step : upper
    % used by multistep methods as a "previous" point:
    sd = feval(model, s, t, param);
    % It will be stored into H
    % First shift the matrix to the right,...
    if (N > 2)
        H = circshift(H, [0, STATE_COLS]);
    end %if
    
    % and overwrite the left part of H with sd1: 
    H(:, 1:STATE_COLS) = sd;
    
    k1 = sd * t_step;
    k2 = feval(model, s+0.5*k1, t+0.5*t_step, param) * t_step;
    k3 = feval(model, s+0.5*k2, t+0.5*t_step, param) * t_step;
    k4 = feval(model, s+k3, t+t_step, param) * t_step;

    s = s + (k1 + 2*k2 + 2*k3 + k4) / 6;

    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
end %for

end %function
