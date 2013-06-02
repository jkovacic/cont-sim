% Solve a set of ordinary differential equations using the
% Hamming's predictor - corrector method. The method is described at:
% http://www-solar.mcs.st-and.ac.uk/~steveb/course/notes/set6_2.pdf
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
%   inputf - name of the function that returns the external input at the specified time
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps


function output = integ_hamming(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

% Dimensions of a matrix of states. Even though it is recommended to be a vertical
% vector (n x 1), the implementation is robust enough to handle other dimensions as well.
[IGNORED, STATE_COLS] = size(initial_condition);

% Check some simulation parameters first
check_sim_params(t_start, t_stop, t_step);

% Check validity of the vector/matrix of states and coefficients:
if ( STATE_COLS <=0 )
    error('Invalid dimensions of initial_condition.');
end %if

upper_limit = t_start + 2 * t_step;
% Note that aux_rk4 will stop at t = t_start + 3*t_step, resulting in
% 4 points, including the initial condition, referred later as (k-3):k

% Make sure that not too many points will be calculated if t_step is too large
if (upper_limit > (t_stop-t_step) )
    upper_limit = t_stop-t_step;
end %if

% The method is not self-starting, so the initial values must be calculated
% using another method. The 4th order Runge - Kutta method is chosen.
[outrk4, S, H] = aux_rk4(model, initial_condition, t_start, upper_limit, t_step, inputf, outputf, param, 4);
s = S(:, 1:STATE_COLS);

% To improve efficiency, preallocate the buffer for output:
[O_ROWS, idx] = size(outrk4);
output = aux_preallocate_output(O_ROWS, t_start, t_stop, t_step);
% Fill the 'output' with initial values
output(:, 1:idx) = outrk4;

% Start of the Hamming's method

ut = feval(inputf, upper_limit+t_step);

% Current index within 'output'
idx = idx+1;
for t = upper_limit+t_step : t_step : t_stop-t_step

    % Predictor (the same as at the Milne - Simpson's method):
    sd = feval(model, s, ut, t, param);
    p = S(:, (3*STATE_COLS+1) : (4*STATE_COLS) ) + 4 * t_step * ...
        ( 2*H(:, (STATE_COLS+1) : (2*STATE_COLS) ) -H(:, 1:STATE_COLS ) + 2*sd ) / 3;
          
    % Corrector:
    ut = feval(inputf, t+t_step);
    pd = feval(model, p, ut, t+t_step, param);
    s = ( 9 * s - S(:, (2*STATE_COLS+1) : (3*STATE_COLS) ) ) / 8 + ...
        3 * t_step * ( pd + 2*sd - H(:, 1:STATE_COLS) ) / 8;
    
    % Shift the history matrices to the right,...
    H = circshift(H, [0, STATE_COLS]);
    S = circshift(S, [0, STATE_COLS]);
    % and overwrite their left parts with sd and s, respectively: 
    H(:, 1:STATE_COLS) = sd;
    S(:, 1:STATE_COLS) = s;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, ut, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    
    % update 'idx'
    idx = idx+1;
end %for

end %function 
