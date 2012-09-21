% A demo script that prepares all necessary data and runs a simulation of a
% pendulum on a movable support as described at:
% http://en.wikipedia.org/wiki/Lagrangian_mechanics#Pendulum_on_a_movable_support

% Vector of model's internal states (must be a 4x1 vector):
% s(1) = x [m]          - position of pendulum base
% s(2) = xd [m/s]       - dx/dt, speed of the pendulum's base
% s(3) = theta [rad]    - pendulum angle
% s(4) = thetad [rad/s] - dtheta/dt, rotational speed of the pendulum

% Vector of desired outputs (a 2x1 vector):
% out(1) = x [cm]
% out(2) = theta [deg]

% Parameters (dimensions) of the model (param can be any vector or matrix):
% param(1) = M [kg]    - mass of the body
% param(2) = m [kg]    - mass of the pendulum
% param(3) = l [m]     - length of the pendulum
% param(4) = g [m/s^2] - gravitational acceleration

% Some sensible values of params:
param = [1, 0.5, 0.5, 9.8067];

% Initial condition at t = 0:
s0 = [0.3, 0, 30*pi/180, 0]';


% Several ODE solving methods are implemented.
% Uncomment the desired one:

%integ_method = "integ_euler";    % Euler method
%integ_method = "integ_ab2";      % 2-step Adams - Bashforth method
%integ_method = "integ_ab3";      % 3-step Adams - Bashforth method
%integ_method = "integ_ab4";      % 4-step Adams - Bashforth method
%integ_method = "integ_ab5";      % 5-step Adams - Bashforth method
integ_method = "integ_rk4";      % 4th order Runge - Kutta method

% Run a simulation cycle t = 0 to 10 s, step = 0.01 s: 
out = feval(integ_method, "model_movable_pendulum", s0, 0, 10, 0.01, "output_movable_pendulum", param);

% Display the results in two plots, one for each output
figure(1);
plot(out(1,:), out(2,:));
xlabel("time  (s)");
ylabel("position  (cm)");

figure(2);
plot(out(1,:), out(3,:));
xlabel("time  (s)");
ylabel("angle  (deg)");
