classdef MilneSimpson < IntegMultiStepAb
    % A differential equation solver, based on the the
    % Milne - Simpson predictor - corrector method. The method is described at:
    % http://math.fullerton.edu/mathews/n2003/milnesimpson/MilneSimpsonProof.pdf
    
    
    methods
        function obj = MilneSimpson(model, input, t_start, t_stop, t_step, init_condition)
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
                error('Invalid dimensions of initial_condition.');
            end %if

            upper_limit = self.t_start + 2 * self.t_step;
            % Note that starter() will stop at t = t_start + 3*t_step, resulting in
            % 4 points, including the initial condition, referred later as (k-3):k

            % Make sure that not too many points will be calculated if t_step is too large
            if (upper_limit > (self.t_stop-self.t_step) )
                upper_limit = self.t_stop-self.t_step;
            end %if

            % The method is not self-starting, so the initial values must be calculated
            % using another method. The 4th order Runge - Kutta method is chosen.
            [outrk4, S, H] = self.starter(upper_limit, 4);
            s = S(:, 1:STATE_COLS);

            % To improve efficiency, preallocate the buffer for output:
            [O_ROWS, idx] = size(outrk4);
            output = self.preallocateOutput(O_ROWS);
            % Fill the 'output' with initial values
            output(:, 1:idx) = outrk4;

            % Start of the Milne - Simpson method:

            ut = self.input.u(upper_limit+self.t_step);

            % Current index within 'output'
            idx = idx+1;
            for t = upper_limit+self.t_step : self.t_step : self.t_stop-self.t_step

                % Predictor:
                sd = self.model.deriv(s, ut, t);
                p = S(:, (3*STATE_COLS+1) : (4*STATE_COLS) ) + 4 * self.t_step * ...
                    ( 2*H(:, (STATE_COLS+1) : (2*STATE_COLS) ) -H(:, 1:STATE_COLS ) + 2*sd ) / 3;

                % Corrector:
                ut = self.input.u(t+self.t_step);
                pd = self.model.deriv(p, ut, t+self.t_step);
                s = S(:, (STATE_COLS+1) : (2*STATE_COLS) ) + self.t_step * ...
                    ( H(:, 1:STATE_COLS ) + 4*sd + pd ) / 3;

                % Shift the history matrices to the right,...
                H = circshift(H, [0, STATE_COLS]);
                S = circshift(S, [0, STATE_COLS]);
                % and overwrite their left parts with sd and s, respectively: 
                H(:, 1:STATE_COLS) = sd;
                S(:, 1:STATE_COLS) = s;

                % Past this point, s represents states at the next point in time, i.e. at t+t_step.
                % This should be kept in mind when calculating output values and applying their time stamp.
                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];

                % update 'idx'
                idx = idx+1;
            end %for           
        end % function
    end % methods
    
end % classdef
