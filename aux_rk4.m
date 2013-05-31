% A convenience function to calculate the first N solutions of a set of ordinary differential
% equations, typically used as a start of multi step methods (e.g. the Adams family).
% This function should not be called directly

% The function is similar to integ_rk4.m. However, this implementation is suited to
% multi step methods, for instance, it also maintains a buffer of the last N values
% of model(sm, tn), etc.

% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   upper_limit - stop time of the simulation run
%   t_step - fixed time step
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%   N - number of previous states to be stored in H
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps
%   S - The last N  states
%   H - buffer of values of model(sn, tn) for the last N points

function [output, S, H] = aux_rk4(model, initial_condition, t_start, upper_limit, t_step, outputf, param, N)

[STATE_ROWS, STATE_COLS] = size(initial_condition);

% A buffer to store previous N-1 values of model(sn, tn), 
% required by multi step methods
H = zeros(STATE_ROWS, (N-1) * STATE_COLS);
S = zeros(STATE_ROWS, N * STATE_COLS);
% Note that some multi step methods (e.g. Milne - Simpson method) also require
% history of states, which is packed into S

% The first set of output values at t = t_start:
initval = feval(outputf, initial_condition, t_start, param);

% To improve efficiency, preallocate the buffer for output:
[O_ROWS, IGNORED] = size(initval);
output = aux_preallocate_output(O_ROWS+1, t_start, upper_limit+t_step, t_step);
% Fill the 'output' with initial values
output(:, 1) = [t_start; initval];

s = initial_condition;
S(:, 1:STATE_COLS) = s;

% Current index within 'output'
idx = 2;
for  t = t_start : t_step : upper_limit
    % used by multistep methods as a "previous" point:
    sd = feval(model, s, t, param);
    % It will be stored into H
    % First shift the matrix to the right,...
    if (N > 2)
        H = circshift(H, [0, STATE_COLS]);
    end %if
    % Also shift the matrix S,  but update it later when s is actually calculated
    if (N > 1)
        S = circshift(S, [0, STATE_COLS]);
    end %if
    
    % and overwrite the left part of H with sd1: 
    H(:, 1:STATE_COLS) = sd;
    
    k1 = sd * t_step;
    k2 = feval(model, s+0.5*k1, t+0.5*t_step, param) * t_step;
    k3 = feval(model, s+0.5*k2, t+0.5*t_step, param) * t_step;
    k4 = feval(model, s+k3, t+t_step, param) * t_step;

    s = s + (k1 + 2*k2 + 2*k3 + k4) / 6;
    % Now the left part of S can be overwritten:
    S(:, 1:STATE_COLS) = s;
    
    val = feval(outputf, s, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    % update 'idx'
    idx = idx+1;
end %for

end %function
