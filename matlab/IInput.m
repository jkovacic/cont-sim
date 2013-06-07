classdef (Abstract) IInput
    % An interface (a pure abstract class with no implemented methods)
    % that all classes for external input must be derived from.
    %
    % A convenience derived class 'InputNan' is provided for cases
    % when no input signal is necessary.
        
    methods
        % An abstract method (must be implemented by derived classes)
        % that returns the value of the input function at the moment 't'
        %
        % Input:
        %   t - moment in time
        % Output:
        %   out - value of the input signal at the moment 't' (may be a vector for multivariable systems)
        out = u(self, t)
    end % properties
    
end % classdef
