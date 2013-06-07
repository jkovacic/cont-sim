classdef InputNan < IInput
    % An implementation of IInput that returns NaN for any 't'.
    % 
    % As such, the class is suitable as a "dummy" class, passed to
    % differential equation solving algorithms when a model does not use
    % any external input.
        
    methods
        function out = u(~, ~)
            % An implementation of u(t) that always returns NaN.
            %
            % Input:
            %   t - moment in time (ignored by this method)
            % Output:
            %   out - NaN
            
            out = NaN;
        end % function
        
    end %methods
    
end %classdef
