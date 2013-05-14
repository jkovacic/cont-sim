% Solve a set of ordinary differential equations using the
% 1st order Nystroem method. The method is described at:
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
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps


function output = integ_nystrom1(model, initial_condition, t_start, t_stop, t_step, outputf, param)

% Dimensions of a matrix of states. Even though it is recommended to be a vertical
% vector (n x 1), the implementation is robust enough to handle other dimensions as well.
[STATE_ROWS, STATE_COLS] = size(initial_condition);

% Check some simulation parameters first
check_sim_params(t_start, t_stop, t_step);

% Check validity of the vector/matrix of states and coefficients:
if ( STATE_COLS <=0 )
    error('Invalid dimensions of initial_condition.');
end %if

upper_limit = t_start;
% Note that aux_rk4 will stop at t = t_start + t_step, resulting in
% 2 points, including the initial condition, referred later as (k-1):k

% Make sure that not too many points will be calculated if t_step is too large
if (upper_limit > (t_stop-t_step) )
    upper_limit = t_stop-t_step;
end %if

% The method is not self-starting, so the initial values must be calculated
% using another method. The 4th order Runge - Kutta mehod is chosen.
[output, S, H] = aux_rk4(model, initial_condition, t_start, upper_limit, t_step, outputf, param, 2);
s = S(:, 1:STATE_COLS);

% Start of the Nystroem's method
%S
for t = upper_limit+t_step : t_step : t_stop-t_step

    sd = feval(model, s, t, param);
    s = S(:, (STATE_COLS+1) : (2*STATE_COLS) ) + 2 * t_step * sd;
          
    % Shift the history matrix to the right,...
    S = circshift(S, [0, STATE_COLS]);
    % and overwrite their left parts with s: 
    S(:, 1:STATE_COLS) = s;
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calcualating output values and applyng their time stamp.
    val = feval(outputf, s, t+t_step, param);
    output = [output, [t+t_step; val]];
end %for

end %function
 
