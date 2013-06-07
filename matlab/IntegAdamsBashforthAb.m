classdef (Abstract) IntegAdamsBashforthAb < IntegMultiStepAb
    % An abstract class with common fuctionality for Adams - Bashforth
    % family of methods for numerical solving of differential equations. 
    %
    % It contains a common engine for all methods of this family.
    %
    % All classes based on Adams - Bashforth methods should be derived from
    % this one.
    
    properties (Access=protected)
        ab_coef;    % A vector of coefficients for the desired Adams - Bashforth method
    end % properties
    
    methods (Access=protected)
        function obj = IntegAdamsBashforthAb(model, input, t_start, t_stop, t_step, init_condition)
            % Constructor.
            %
            % Input:
            %   model - model of the simulated system (must be an implementation of IModel)
            %   input - external input(s) to the system (must be an implementation of IInput),
            %           if it is not necessary, InputNan may be passed
            %   t_start - start time of the simulation run
            %   t_stop - stop time of the simulation run
            %   t_step - fixed time step
            %   init_condition - vector of model's states' values at t=t_start
            
            obj = obj@IntegMultiStepAb(model, input, t_start, t_stop, t_step, init_condition);
        end % function
    end % methods
        
    methods
        function output = run(self)
            % Start the simulation run.
            % The common "engine" is used for all Adams - Bashforth
            % methods.
            %
            % Output:
            %   output - matrix of output values (as defined by model.out), prepended by time stamps
            
            % number of steps:
            N = numel(self.ab_coef);

            % Dimensions of a matrix of states. Even though it is recommended to be a vertical
            % vector (n x 1), the implementation is robust enough to handle other dimensions as well.
            [~, STATE_COLS] = size(self.initial);

            % Check some simulation parameters first
            self.checkSimParams();

            % Check validity of the vector/matrix of states and coefficients:
            if ( STATE_COLS <=0 )
                error('Invalid dimensions of initial_condition.');
            end %if

            if ( N < 2 )
                error('Only 2- or more- steps Adams - Bashforth methods are supported.');
            end %if


            upper_limit = self.t_start + (N-1)*self.t_step;

            % Make sure that not too many points will be calculated if t_step is too large
            if (upper_limit > (self.t_stop-self.t_step) )
                upper_limit = self.t_stop-self.t_step;
            end %if

            % The method is not self-starting, so the initial values must be calculated
            % using another method. The 4th order Runge - Kutta method is chosen.
            [outrk4, S, H] = self.starter(upper_limit, N);
            s = S(:, 1:STATE_COLS);

            % To improve efficiency, preallocate the buffer for output:
            [O_ROWS, idx] = size(outrk4);
            output = self.preallocateOutput(O_ROWS);

            % Fill the 'output' with initial values
            output(:, 1:idx) = outrk4;

            % Now the Adams - Bashforth method can start

            ut = self.input.u(upper_limit+self.t_step);

            % Current index within 'output'
            idx = idx+1;
            for t = upper_limit+self.t_step : self.t_step : self.t_stop-self.t_step,
                sd = self.model.deriv(s, ut, t);
                s = s + self.t_step * sd * self.ab_coef(1);
                for i = 2:N
                    s = s + self.t_step * self.ab_coef(i) * ...
                        H(:, ((i-2)*STATE_COLS+1) : ((i-1)*STATE_COLS) );
                end %for

                % Shift the H to the right,...
                if (N > 2)
                    H = circshift(H, [0, STATE_COLS]);
                end %if
                % and overwrite the left part of H with sd: 
                H(:, 1 : STATE_COLS) = sd;

                % Past this point, s represents states at the next point in time, i.e. at t+t_step.
                % This should be kept in mind when calculating output values and applying their time stamp.
                ut = self.input.u(t+self.t_step);
                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];

                % update 'idx'
                idx = idx+1;
            end % for            
        end % function
        
    end % methods
    
end % classdef
