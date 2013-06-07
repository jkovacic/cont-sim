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

% Instantiate the class with implementation of the model:
m = MovablePendulum(param);

% Since no external input is necessary, a "dummy" class (that returns NaN)
% can be passed to the selected diff. equation solving algorithms
inp = InputNan();

% States' initial condition at t = 0:
s0 = [0.3, 0, 30*pi/180, 0]';


% Parameters of a simulation cycle: t = 0 to 10 s, step = 0.01 s:
tst = 0;
tstop = 10;
step = 0.01;

% Several ODE solving methods are implemented.
% Uncomment the desired one:

%integ_method = 'Euler';                % Euler method
%integ_method = 'Nystrom1';             % 1st order Nystroem's method
%integ_method = 'AB2';                  % 2-step Adams - Bashforth method
%integ_method = 'AB3';                  % 3-step Adams - Bashforth method
%integ_method = 'AB4';                  % 4-step Adams - Bashforth method
%integ_method = 'AB5';                  % 5-step Adams - Bashforth method
%integ_method = 'ABM1';                 % 1-step Adams - Bashforth - Moulton predictor - corrector method
%integ_method = 'ABM2';                 % 2-step Adams - Bashforth - Moulton predictor - corrector method
%integ_method = 'ABM3';                 % 3-step Adams - Bashforth - Moulton predictor - corrector method
%integ_method = 'ABM4';                 % 4-step Adams - Bashforth - Moulton predictor - corrector method
%integ_method = 'ABM5';                 % 5-step Adams - Bashforth - Moulton predictor - corrector method
%integ_method = 'MilneSimpson';         % Milne - Simpson predictor - corrector method
%integ_method = 'MilneSimpsonImproved'; % improved Milne - Simpson predictor - corrector method
%integ_method = 'Hamming';              % Hamming's predictor - corrector method 
%integ_method = 'HammingImproved';      % improved Hamming's predictor - corrector method 
%integ_method = 'RK2';                  % 2nd order Runge - Kutta method (alpha = 1/2, the midpoint method)
%integ_method = 'Heun2';                % 2nd order Heun method (i.e. 2nd order Runge - Kutta, where alpha=1)
%integ_method = 'Ralston2';             % 2nd order Ralston's method (i.e. 2nd order Runge - Kutta, where alpha=2/3)
%integ_method = 'RK3';                  % 3rd order Runge - Kutta method
%integ_method = 'Heun3';                % 3rd order Heun method
integ_method = 'RK4';                  % 4th order Runge - Kutta method
%integ_method = 'RK_3_8';               % Runge - Kutta 3/8 method
%integ_method = 'Ralston4'              % 4th order Ralston's method
%integ_method = 'Gill4';                % 4th order Gill method
%integ_method = 'Nystrom5';             % 5th order Nystroem's method
%integ_method = 'Butcher6';             % 6th order Butcher's method
%integ_method = 'Verner8';              % 8th order Verner's method


% Instantiate the diff. equation solving class:
eval(strcat('eng = ', integ_method, '(m, inp, tst, tstop, step, s0);'));
% it should evaluate into something like this:
% eng = RK4(m, inp, tst, tstop, step, s0);

% And finally run a simulation cycle:
out = eng.run();

% Display the results in two plots, one for each output
figure(1);
plot(out(1,:), out(2,:));
xlabel('time (s)');
ylabel('position (cm)');

figure(2);
plot(out(1,:), out(3,:));
xlabel('time (s)');
ylabel('angle (deg)');
