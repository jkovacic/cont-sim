% Elements of a general Butcher tableau for 2nd order Runge - Kutta methods.
% For more details, see:
% http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Second-order_methods_with_two_stages
%
% The function should not be called directly

% Input:
%   alpha - parameter of a 2nd order Runge - Kutta methods
%
% Output:
%   A - matrix A of the Butcher tableau
%   B - vector B of the Butcher tableau
%   C - vector C of the Butcher tableau

function [A, B, C] = aux_rk2_butcher(alpha)

if ( abs(alpha) <= eps )
    error('Alpha cannot be zero!');
end %if

A = zeros(2, 2);
B = zeros(1, 2);
C = zeros(2, 1);

% for a slight improvement in speed, calculate 1/(2*alpha) only once:
invalpha = 1/(2*alpha);

A(2, 1) = alpha;
C(2) = alpha;

B(1) = 1 - invalpha;
B(2) = invalpha;

end %function
