% A function to calculate desired outputs from the states space representation 
% of a pendulum on a movable support as derived at:
% http://en.wikipedia.org/wiki/Lagrangian_mechanics#Pendulum_on_a_movable_support
% 
% The function will typically be executed by differential equation solving functions,
% so its parameters must conform to the model "interface" as described in README.
%
% Input:
%   s - vector of values of states at the specified moment in time
%   t - moment in time (not used by this function)
%   param - values of model parameters (not used by this function)
%
% Output:
%   out - vector of desired output values at the specified time
%
% The desired output vector consists of the following values:
% x [cm]      - horizontal position of the pendulum support 
% theta [deg] - inclination of the pendulum (relative to the vertical line)
%
% The model's states have the following meanings:
% s(1) = x [m]          - horizontal position of the pendulum support 
% s(2) = xp [m/s]       - dx/dt, horizontal speed of the pendulum support
% s(3) = theta [rad]    - inclination of the pendulum (relative to the vertical line)
% s(4) = thetap [rad/s] - dtheta/dt, rotational speed of the pendulum

function out = output_movable_pendulum(s, t, param)


% Vector of desired outputs ( a 2x1 vector):
% out(1) = x [cm]
% out(2) = theta [deg]
%
% The output is derived from the following vector of states s:
% s(1) = x [m]          - position of pendulum base
% s(2) = xp [m/s]       - dx/dt, speed of the pendulum's base
% s(3) = theta [rad]    - pendulum angle
% s(4) = thetap [rad/s] - dtheta/dt, rotational speed of the pendulum


N_STATES = 4;

if ( N_STATES ~= numel(s) )
    error("Invalid number of states, expected %d", N_STATES);
end %if

% The expected number of parameters is also 4, however in no relation to N_STATES
if ( numel(param) < 4 )
    error("Invalid number of model's parameters");
end % if


% Note that x must be coverted from meters to centimeters (multiplying by 100)
% and theta must be coverted from radians to degrees (multiplying by 180/pi).
out = zeros(2,1);
out(1) = s(1) * 100;
out(2) = s(3) * 180/pi;

% Could also be impleneted as:
% out = [s(1)*100, s(3)*180/pi]';

end %function
