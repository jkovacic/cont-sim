% Solve a set of ordinary differential equations using the 3rd order Heun method 
% (similar to the 3rd order Runge - Kutta method). The method is described at:
% http://de.wikipedia.org/wiki/Runge-Kutta-Verfahren#Beispiele
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

function output = integ_heun3(model, initial_condition, t_start, t_stop, t_step, outputf, param)

% Elements of the Butcher tableau:
A = [0, 0, 0; 1, 0, 0; 0, 2, 0] / 3;
B = [1, 0, 3] / 4;
C = [0, 1, 2]' / 3;

% passed to the general implementation of explicit Runge - Kutta methods
output = aux_rk_expl(model, initial_condition, t_start, t_stop, t_step, outputf, param, A, B, C);

end % function  
 
