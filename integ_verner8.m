% Solve a set of ordinary differential equations using the 8th order Verner's method.
% The method is described at:
% http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_verner.c
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

function output = integ_verner8(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param)

s21 = sqrt(21);

% Elements of the Butcher tableau:
A = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
      0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
      0.25, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
      1/7, (-7-3*s21)/98, (21+5*s21)/49, 0, 0, 0, 0, 0, 0, 0, 0; ...
      (11+s21)/84, 0, (18+4*s21)/63, (21-s21)/252, 0, 0, 0, 0, 0, 0, 0; ...
      (5+s21)/48, 0, (9+s21)/36, (-231+14*s21)/360, (63-7*s21)/80, 0, 0, 0, 0, 0, 0; ...
      (10-s21)/42, 0, (-432+92*s21)/315, (633-145*s21)/90, (-504+115*s21)/70, (63-13*s21)/35, 0, 0, 0, 0, 0; ...
      1/14, 0, 0, 0, (14-3*s21)/126, (13-3*s21)/63, 1/9, 0, 0, 0, 0; ...
      1/32, 0, 0, 0, (91-21*s21)/576, 11/72, (-385-75*s21)/1152, (63+13*s21)/128, 0, 0, 0; ...
      1/14, 0, 0, 0, 1/9, (-733-147*s21)/2205, (515+111*s21)/504, (-51-11*s21)/56, (132+28*s21)/245, 0, 0; ...
      0, 0, 0, 0, (-42+7*s21)/18, (-18+28*s21)/45, (-273-53*s21)/72, (301+53*s21)/72, (28-28*s21)/45, (49-7*s21)/18, 0 ];

B = [9, 0, 0, 0, 0, 0, 0, 49, 64, 49, 9] / 180;
C = [ 0, 0.5, 0.5, (7+s21)/14, (7+s21)/14, 0.5, ...
      (7-s21)/14, (7-s21)/14, 0.5, (7+s21)/14, 1]';

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, inputf, outputf, param, A, B, C);

end % function  
