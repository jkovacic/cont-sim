% Checks validity of simulation cycle parameters,
% e.g., all parameters must be scalars, t_stop must be greater than t_start,
% t_step must be positive, etc. If any parameter is invalid, a short
% description of the problem will be displayed and an error will be thrown.
%
% t_start - desired start time of the simulation run
% t_stop - desired stop time of the simulation cycle
% t_step - desired time step

function check_sim_params(t_start, t_stop, t_step)


if ( numel(t_start) ~= 1 )
    error ("Start time must be a scalar");
end %if

if ( numel(t_stop) ~= 1 )
    error ("Stop time must be a scalar");
end %if

if ( numel(t_step) ~= 1 )
    error ("Time step must be a scalar");
end %if

if ( t_stop < t_start )
    error ("Stop time must be greater than start time");
end %if

if ( t_step <= 0 )
    error("t_step must bepositive");
end % if

end %function
