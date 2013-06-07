% Preallocates the buffer for output of fixed step algorithms
% The function should be called internally only.

% Input:
%   rows - number of output's rows
%   t_start - start time of the simulation run
%   t_stop - stop time of the simulation run
%   t_step - fixed time step  
%
% Output:
%   output - a [rows, floor((t_stop-t_start)/t_step)+1] matrix, filled with Nan 

function output = aux_preallocate_output(rows, t_start, t_stop, t_step)
    % It is assumed that regularity of parameters is checked 
    % by the caller beforehand
    
    output = NaN(rows, floor((t_stop-t_start)/t_step)+1);
end %function
