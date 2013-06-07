% A constant function that always returns NaN.
% As such, the function is suitable as a "dummy" function, passed to
% differential equation solving algorithms when a model does not use
% any external input.
%
% The function will typically be executed by differential equation solving functions,
% so its parameters must conform to the model "interface" as described at
% https://github.com/jkovacic/cont-sim/wiki/Basic-instructions.
%
% Input:
%   t - moment in time (not used by this function)
%
% Output:
%   out - functions value at the moment of 't' (NaN in this function)

function out = input_nan(t)
    out = NaN;
end % function
