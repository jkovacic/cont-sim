% A general function to solve a set of ordinary differential equations using
% one of Adams - Bashforth methods.
% This function should only be called from integ_ab?.m

% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step
%   inputf - name of the function that returns the external input at the specified time
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%   ab_coefficients - a vector of coefficients for the desired Adams - Bashforth method
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps

function output = aux_ab_general(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, ab_coefficients) 

% number of steps:
N = numel(ab_coefficients);

% Dimensions of a matrix of states. Even though it is recommended to be a vertical
% vector (n x 1), the implementation is robust enough to handle other dimensions as well.
[IGNORED, STATE_COLS] = size(initial_condition);


% Check some simulation parameters first
check_sim_params(t_start, t_stop, t_step);

% Check validity of the vector/matrix of states and coefficients:
if ( STATE_COLS <=0 )
    error('Invalid dimensions of initial_condition.');
end %if

if ( N < 2 )
    error('Only 2- or more- steps Adams - Bashforth methods are supported.');
end %if


upper_limit = t_start + (N-1)*t_step;

% Make sure that not too many points will be calculated if t_step is too large
if (upper_limit > (t_stop-t_step) )
    upper_limit = t_stop-t_step;
end %if

% The method is not self-starting, so the initial values must be calculated
% using another method. The 4th order Runge - Kutta method is chosen.
[outrk4, S, H] = aux_rk4(model, initial_condition, t_start, upper_limit, t_step, inputf, outputf, param, N);
s = S(:, 1:STATE_COLS);

% To improve efficiency, preallocate the buffer for output:
[O_ROWS, idx] = size(outrk4);
output = aux_preallocate_output(O_ROWS, t_start, t_stop, t_step);

% Fill the 'output' with initial values
output(:, 1:idx) = outrk4;

% Now the Adams - Bashforth method can start

ut = feval(inputf, upper_limit+t_step);

% Current index within 'output'
idx = idx+1;
for t = upper_limit+t_step : t_step : t_stop-t_step,
    sd = feval(model, s, ut, t, param);
    s = s + t_step * sd * ab_coefficients(1);
    for i = 2:N
        s = s + t_step * ab_coefficients(i) * H(:, ((i-2)*STATE_COLS+1) : ((i-1)*STATE_COLS) );
    end %for
    
    % Shift the H to the right,...
    if (N > 2)
        H = circshift(H, [0, STATE_COLS]);
    end %if
    % and overwrite the left part of H with sd: 
    H(:, 1 : STATE_COLS) = sd;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    ut = feval(inputf, t+t_step);
    val = feval(outputf, s, ut, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    
    % update 'idx'
    idx = idx+1;
end % for

end % function
