% State space representation of the model of a pendulum on a movable support as derived at:
% http://en.wikipedia.org/wiki/Lagrangian_mechanics#Pendulum_on_a_movable_support
%
% The function will typically be executed by differential equation solving functions,
% so its parameters must conform to the model "interface" as described at
% https://github.com/jkovacic/cont-sim/wiki/Basic-instructions.
%
% Input:
%   s - vector of values of states at the specified moment in time
%   t - moment in time (not used by this function)
%   param - values of model parameters
%
% Output:
%   sd - vector of derivatives of states with respect of time at the time = t
%
% Instead of being hardcoded in this file, the model's parameters (dimensions) are passed
% in param, allowing to test responses of the model using different parameter values.
% The model expects the parameters in a vector as defined below:
% 
% param(1) = M [kg]    - mass of the movable body
% param(2) = m [kg]    - mass of the pendulum point
% param(3) = l [m]     - length of the pendulum
% param(4) = g [m/s^2] - gravitational acceleration
%
% The model's states have the following meanings:
% s(1) = x [m]          - horizontal position of the pendulum support 
% s(2) = xd [m/s]       - dx/dt, horizontal speed of the pendulum support
% s(3) = theta [rad]    - inclination of the pendulum (relative to the vertical line)
% s(4) = thetad [rad/s] - dtheta/dt, rotational speed of the pendulum


function sd = model_movable_pendulum(s, t, param)


% Required number of states to check the correctness of given 's':
N_STATES = 4;

if ( N_STATES ~= numel(s) )
    error('Invalid number of states, expected %d', N_STATES);
end %if

% The expected number of parameters is also 4, however in no relation to N_STATES
if ( numel(param) < 4 )
    error('Invalid number of model''s parameters');
end % if

% Elements of params are assigned to variables with meaningful names.
% This consumes a little bit more memory, however it siginificantly improves readability and
% maintainability of the code and reduces chances of errors 
% (e.g. accidentally passing a wrong index to param)
m = param(1);
M = param(2);
l = param(3);
g = param(4);

% For exactly the same reasons, values of states are assigned to variables with meaningful names: 
x = s(1);        % [m]     - horizontal position of the pendulum support 
xd = s(2);       % [m/s]   - dx/dt, horizontal speed of the pendulum support
theta = s(3);    % [rad]   - inclination of the pendulum (relative to the vertical line)
thetad = s(4);   % [rad/s] - dtheta/dt, rotational speed of the pendulum


% The derived model is a set of two ordinary differential equations:
%
% (M+m)*xdd + m*l*cos(theta)*thetadd - m*l*sin(theta)*thetad^2 = 0
% thetadd + xdd*cos(theta)/l + g*sin(theta)/l = 0                      (1)
%
% The highest order derivatives of both states (xdd and thetadd) must be expressed
% as explicite functions of the lower level derivatives, inputs, time, params, etc.
%
% From the first equation, xdd can be derived as:
% xdd = m*l*(thetad^2*sin(theta)-thetadd*cos(theta))/(m+M)             (2)
%
% and put into the the second equation which turns into:
% thetadd + m*(thetad^2*sin(theta)*cos(theta)-thetadd*cos^2(theta))/(m+M) +g*sin(theta)/l = 0       (3)
%
% thetadd can be derived from it as:
% thetadd = (-sin(theta)*(m*thetad^2*cos(theta)/(m+M))+g/l)/(1-m*cos^2(theta)/(m+M))                (4)
%
% Put the thetadd into (2) and xdd can be expressed as:
% xdd = (m*l/(m+M))*(thetad^2*sin(theta)+(sin(theta)*cos(theta)*(m*thetad^2*cos(theta)/(m+M))+g/l)/ \\
%   (1-m*cos^2(theta)/(m+M)))                                                                       (5)
%
% The equations (4) and (5) can be used to calculate sd.
% Note that (1-m*cos^2(theta)/(m+M)) can never limit to zero unless M is negligible comparing to m.
% As this is unlikely in typical situations, the division by zero is very unlikely
% and will not be handled separately


% Sine and cosine of theta are used often in the model. Since trigonometrical functions are 
% computationally complex, store the theta's sine and cosine into temporary values:
sth = sin(theta);
cth = cos(theta);

% m/(m+M) and g/l also appear often, so the values will be stored into separate variables:
mrel = m/(m+M);
gl = g/l;

% Finally calculate numerical values of both second order derivatives as derived above:
xdd = mrel*l*(thetad*thetad*sth + (sth*cth*(mrel*thetad*thetad*cth+gl)) / (1-mrel*cth*cth));
thetadd = -sth*(mrel*thetad*thetad*cth+gl)/(1-mrel*cth*cth);

% And appropriately fill the output vector of derivatives:
sd = s;    % to preserve the dimensions and avoid errors due to addition of incompatible vectors/matrices
sd(1) = xd;
sd(2) = xdd;
sd(3) = thetad;
sd(4) = thetadd;

end %function
