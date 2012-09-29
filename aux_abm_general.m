% A general function to solve a set of ordinary differential equations using
% one of Adams - Bashforth methods.
% This function should only be called from integ_abm?.m

% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%   ab_coefficients - a vector of coefficients for the desired Adams - Bashforth predictor method
%   am_coefficients - a vector of coefficients for the desired Adams - Moulton corrector method
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps

function output = aux_abm_general(model, initial_condition, t_start, t_stop, t_step, outputf, param, ab_coefficients, am_coefficients) 

% number of steps of the Adams - Bashforth method:
N = numel(ab_coefficients);
% number of steps of the Adams - Moulton method:
Nam = numel(am_coefficients);

% Dimensions of a matrix of states. Even though it is recommended to be a vertical
% vector (n x 1), the implementation is robust enough to handle other dimensions as well.
[STATE_ROWS, STATE_COLS] = size(initial_condition);


% Check some simulation parameters first
check_sim_params(t_start, t_stop, t_step);

% Check validity of the vector/matrix of states and coefficients:
if ( STATE_COLS <=0 )
    error("Invalid dimensions of initial_condition.");
end %if

if ( N < 2 )
    error("Only 2- or more- steps Adams - Bashforth methods are supported.");
end %if

if ( Nam<1 || Nam >(N+1) )
    error("Invalid number of Adams - Moulton coefficients");
end %if

upper_limit = t_start + (N-1)*t_step;

% Make sure that not too many points will be calculated if t_step is too large
if (upper_limit > (t_stop-t_step) )
    upper_limit = t_stop-t_step;
end %if

% The method is not self-starting, so the initial values must be calculated
% using another method. The 4th order Runge - Kutta mehod is chosen.
[output, S, H] = aux_rk4(model, initial_condition, t_start, upper_limit, t_step, outputf, param, N);
s = S(:, 1:STATE_COLS);

% Now the Adams - Bashforth method can start

for t = upper_limit+t_step : t_step : t_stop-t_step,
    % Predictor method:
    sd = feval(model, s, t, param);
    p = s + t_step * sd * ab_coefficients(1);
    for i = 2:N
        p = p + t_step * ab_coefficients(i) * H(:, ((i-2)*STATE_COLS+1) : ((i-1)*STATE_COLS) );
    end %for
    
    % Adams - Moulton corrector method:
    pd = feval(model, p, t+t_step, param);
    s = s + t_step * pd * am_coefficients(1);
    if ( Nam >=2 )
        s = s + t_step * sd * am_coefficients(2);
    end %if
    for i = 3:Nam
        s = s + t_step * am_coefficients(i) * H(:, ((i-3)*STATE_COLS+1) : ((i-2)*STATE_COLS) );
    end %for
    
    % Shift the H to the right,...
    if (N > 2)
        H = circshift(H, [0, STATE_COLS]);
    end %if
    % and overwrite the left part of H with sd: 
    H(:, 1 : STATE_COLS) = sd;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
    
end % for

end % function 
