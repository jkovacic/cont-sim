classdef Nystrom1 < IntegMultiStepAb
    % A differential equation solver, based on the 1st order Nystroem
    % method. The method is described at:
    % http://www-solar.mcs.st-and.ac.uk/~steveb/course/notes/set6_2.pdf
    
    methods
        function obj = Nystrom1(model, input, t_start, t_stop, t_step, init_condition)
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
        
        function output = run(self)
            % Start the simulation run.
            %
            % Output:
            %   output - matrix of output values (as defined by model.out), prepended by time stamps
            
            % Dimensions of a matrix of states. Even though it is recommended to be a vertical
            % vector (n x 1), the implementation is robust enough to handle other dimensions as well.
            [~, STATE_COLS] = size(self.initial);

            % Check some simulation parameters first
            self.checkSimParams();

            % Check validity of the vector/matrix of states and coefficients:
            if ( STATE_COLS <=0 )
                error('Invalid dimensions of initial condition.');
            end %if

            upper_limit = self.t_start;
            % Note that starter() will stop at t = t_start + t_step, resulting in
            % 2 points, including the initial condition, referred later as (k-1):k

            % Make sure that not too many points will be calculated if t_step is too large
            if (upper_limit > (self.t_stop-self.t_step) )
                upper_limit = self.t_stop-self.t_step;
            end %if

            % The method is not self-starting, so the initial values must be calculated
            % using another method. The 4th order Runge - Kutta method is chosen.
            [outrk4, S, ~] = self.starter(upper_limit, 2); 
            s = S(:, 1:STATE_COLS);

            % To improve efficiency, preallocate the buffer for output:
            [O_ROWS, idx] = size(outrk4);
            output = self.preallocateOutput(O_ROWS);
            % Fill the 'output' with the initial value:
            output(:, 1:idx) = outrk4;

            % Start of the Nystroem's method:

            ut = self.input.u(upper_limit+self.t_step);

            % Current index within 'output':
            idx = idx+1;
            for t = upper_limit+self.t_step : self.t_step : self.t_stop-self.t_step

                sd = self.model.deriv(s, ut, t);
                s = S(:, (STATE_COLS+1) : (2*STATE_COLS) ) + 2 * self.t_step * sd;

                % Shift the history matrix to the right,...
                S = circshift(S, [0, STATE_COLS]);
                % and overwrite their left parts with s: 
                S(:, 1:STATE_COLS) = s;

                % Past this point, s represents states at the next point in time, i.e. at t+t_step.
                % This should be kept in mind when calculating output values and appliyng their time stamp.
                ut = self.input.u(t+self.t_step);
                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];

                % update 'idx'
                idx = idx+1;
            end %for
        end % function
    end % methods
    
end % classdef
