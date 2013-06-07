% Solve a set of ordinary differential equations using the 4th order Ralston's method.
% The method is described at:
% http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_ralston_4.c
%
% Function's input and output parameters should conform to the general integ
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

function output = integ_ralston4(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)
  
s5 = sqrt(5);

% Elements of the Butcher tableau:
A = [0, 0, 0, 0; 0.4, 0, 0, 0; (-2889+1428*s5)/1024, (3785-1620*s5)/1024, 0, 0; ...
     (-3365+2094*s5)/6040, (-975-3046*s5)/2552, (467040+203968*s5)/240845, 0];
B = [(263+24*s5)/1812, (125-1000*s5)/3828, 1024*(3346+1623*s5)/5924787, (30-4*s5)/123];
C = [0, 0.4, 0.875-0.1875*s5, 1]';

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C);

end % function  
