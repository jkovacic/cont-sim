classdef (Abstract) IModel
    % An interface (a pure abstract class with no implemented methods)
    % that all model classes must be derived from.
    %
    % Instance of classes that implement this interface are passed to
    % classes that solve ordinary differential equations and call
    % implementations of abstract methods depending on an algorithm.
   
    properties (Abstract, Constant)
        % A derived model must set this constant to the actual number of
        % all state variables.
        %
        % The qualifier 'Constant' also automatically implies 'Static'.
        N_STATES;
    end % properties
    
    methods
        % An abstract method (must be implemented by derived classes)
        % that returns derivatives of all state variables.
        %
        % Input:
        %   s - a vector of state values at the moment 't' (required)
        %   u - external input at the moment 't' (may be ignored)
        %   t - moment in time (may be ignored)
        % Output:
        %   sd - vector of derivations of state variables
        sd = deriv(self, s, u, t)
        
        % An abstract method (must be implented by derived classes)
        % that returns the model's output signals.
        %
        % Input:
        %   s - vector of state variables at the monet of 't' (required)
        %   u - external input at the moment 't' (may be ignored)
        %   t - moment in time (may be ignored)
        % Output:
        %   output - vector of model's output signals at 't'
        output = out(self, s, u, t)
    end % methods
         
end % classdef
