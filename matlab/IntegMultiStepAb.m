classdef (Abstract) IntegMultiStepAb < IntegAb
    % An abstract class with common functionality for multi step methods,
    % i.e. methods that use values from previous steps, e.g. Adams -
    % Basforth and Adams - Bashforth - Moulton methods, Hamming methods,
    % Milne - Simpson methods, etc.
    %
    % As multi step methods are not self starting, the class contains a
    % method to calculate first steps using a single step method (the 4th order 
    % Runge - Kutta method).
    %
    % All classes based on multi step algorithms should be derived from
    % this one.
    
    methods (Access=protected)
        function obj = IntegMultiStepAb(model, input, t_start, t_stop, t_step, init_condition)
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
            
            obj = obj@IntegAb(model, input, t_start, t_stop, t_step, init_condition);
        end % function
        
        function [output, S, H] = starter(self, upper_limit, N)
            % Calculates the first N solutions of a set of ordinary differential
            % equations, typically used as a start of multi step methods.
            % The 4th order Runge - Kutta method is chosen.
            %
            % Input:
            %   upper_limit - stop time of the simulation run of this method
            %   N - number of previous states to be stored in H
            %
            % Output:
            %   output - matrix of output values, prepended by time stamps
            %   S - The last N  states
            %   H - buffer of values of model(sn, tn) for the last N points
            
            [STATE_ROWS, STATE_COLS] = size(self.initial);

            % A buffer to store previous N-1 values of model(sn, tn), 
            % required by multi step methods
            H = zeros(STATE_ROWS, (N-1) * STATE_COLS);
            S = zeros(STATE_ROWS, N * STATE_COLS);
            % Note that some multi step methods (e.g. Milne - Simpson method) also require
            % history of states, which is packed into S

            % The first set of output values at t = t_start:
            ut = self.input.u(self.t_start);
            initval = self.model.out(self.initial, ut, self.t_start);

            % To improve efficiency, preallocate the buffer for output:
            [O_ROWS, ~] = size(initval);
            output = NaN(O_ROWS+1, floor((upper_limit+self.t_step-self.t_start)/self.t_step)+1);
            
            % Fill the 'output' with initial values
            output(:, 1) = [self.t_start; initval];

            s = self.initial;
            S(:, 1:STATE_COLS) = s;

            % Current index within 'output'
            idx = 2;
            for  t = self.t_start : self.t_step : upper_limit
                % used by multistep methods as a "previous" point:
                sd = self.model.deriv(s, ut, t);
                % It will be stored into H
                % First shift the matrix to the right,...
                if (N > 2)
                    H = circshift(H, [0, STATE_COLS]);
                end %if
                % Also shift the matrix S,  but update it later when s is actually calculated
                if (N > 1)
                    S = circshift(S, [0, STATE_COLS]);
                end %if

                % and overwrite the left part of H with sd1: 
                H(:, 1:STATE_COLS) = sd;

                k1 = sd * self.t_step;
                ut = self.input.u(t+0.5*self.t_step);
                k2 = self.model.deriv(s+0.5*k1, ut, t+0.5*self.t_step) * self.t_step;
                k3 = self.model.deriv(s+0.5*k2, ut, t+0.5*self.t_step) * self.t_step;
                ut = self.input.u(t+self.t_step);
                k4 = self.model.deriv(s+k3, ut, t+self.t_step) * self.t_step;

                s = s + (k1 + 2*k2 + 2*k3 + k4) / 6;
                % Now the left part of S can be overwritten:
                S(:, 1:STATE_COLS) = s;

                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];
                % update 'idx'
                idx = idx+1;
            end %for
        end % function
    end % methods
    
end % classdef
