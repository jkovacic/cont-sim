% A general function to solve a set of ordinary differential equations using
% one of explicit Runge - Kutta methods.
% This function should only be called from integ_rk?.m

% Input:
%   model - name of a function that implements the mathematical model as a state-space representation
%   initial_condition - states at t = t_start
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step
%   inputf - name of the function that returns the external input at the specified time
%   outputf - name of the function that calculates the desired output values from the internal states
%   param - vector parameter values, passed to 'model' and 'outputf'
%   A - matrix A of the Butcher tableau (its diagonal and upper triangle will never be read)
%   B - vector B of the Butcher tableau
%   C - vector C of the Butcher tableau (the first element will automatically be set to 0)
%
% Output:
%   output - vector of output values (as defined by 'outputf'), prepended by time stamps

% More details about elements of a Butcher tableau can befound at:
% http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods

function output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C)

% Check validity of Butcher tableau elements:
% A must be a square matrix:
[N, columns] = size(A);
if ( N ~= columns )
    error('A must be square matrix');
end %if

% B should be a 1xN vector, however any matrix with N elemnts will do
if ( numel(B) ~= N )
    error('Invalid size of B');
end %if

% C should be a Nx1 vector, however any matrix with N elemnts will do
if ( numel(C) ~= N )
    error('Invalid size of C');
end %if

% Some elements of C and A should be 0.
% Instead of checking their values and reporting an error,
% those elements will be set to zeros:
C(1) = 0;
for i = 1:N
    A(i, i:N) = zeros(1, N-i+1);
end %for

% Consistency of the Butcher tableau:
% Note, due to limited accuracy of floating point arithmetics,
% a bit larger threshold for "equality" is necessary for some methods, e.g. Ralston's or Verner's.
crit = 7 * eps;
for i = 1:N
    % This might be useful to determine the threshold
    % diff = sum(A(i, :)) - C(i);
    % printf("diff: %f    crit: %f\n", diff/eps, crit/eps);
    if ( abs( sum(A(i, :)) - C(i) ) > crit )
        error('Inconsistent Butcher tableau');
    end %if
end %for

% Check some simulation parameters
check_sim_params(t_start, t_stop, t_step);

% Dimensions of a matrix of states. Even though it is recommended to be a vertical
% vector (n x 1), the implementation is robust enough to handle other dimensions as well.
[STATE_ROWS, STATE_COLS] = size(initial_condition);

% Check validity of the vector/matrix of states and coefficients:
if ( STATE_COLS <=0 )
    error('Invalid dimensions of initial_condition.');
end %if

% Matrix to store values of kn
K = zeros(STATE_ROWS, N*STATE_COLS);

% The first set of output values at t = t_start:
ut = feval(inputf, t_start);
initval = feval(outputf, initial_condition, ut, t_start, param);

% To improve efficiency, preallocate the buffer for output:
[NOUT, IGNORED] = size(initval);
output = aux_preallocate_output(1+NOUT, t_start, t_stop, t_step);
% Fill the 'output' with the initial value:
output(:, 1) = [t_start; initval];


% Start of the Runge - Kutta algorithm:

s = initial_condition;
% Current index within 'output'
idx = 2;
for t = t_start : t_step : t_stop-t_step
    % The current value of s is used several times inside the for i loop.
    % However, to eleminate the need for another foor loop (after for i),
    % s can be updated inside the same for-loop. As this would break the algorithm
    % for calculation of kn, stemp was introduced.
    stemp = s;
    for i = 1:N
        rktemp = stemp;
        for j = 1: (i-1)
            rktemp = rktemp + A(i, j) * K(:, ((j-1)*STATE_COLS+1) : (j*STATE_COLS) );
        end %for
        
        % ki and its position in K:
        ui = feval(inputf, t+C(i)*t_step);
        K(:, ((i-1)*STATE_COLS+1) : (i*STATE_COLS) ) = ...
            t_step * feval(model, rktemp, ui, t+C(i)*t_step, param);
            
        s = s + B(i) * K(:, ((i-1)*STATE_COLS+1) : (i*STATE_COLS) );
    end %for
    
    % Past this point, s represents states at the next point in time, i.e. at t+t_step.
    % This should be kept in mind when calculating output values and applying their time stamp.
    ut = feval(inputf, t+t_step);
    val = feval(outputf, s, ut, t+t_step, param);
    output(:, idx) = [t+t_step; val];
    % update 'idx'
    idx = idx + 1;
end % for

end %function
