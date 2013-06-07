classdef Euler < IntegAb
    % A differential equation solver, based on the the most simple
    % numerical method, i.e. the Euler method. The method is described at:
    % http://en.wikipedia.org/wiki/Euler_method
    
    methods
        function obj = Euler(model, input, t_start, t_stop, t_step, init_condition)
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
        
        function output = run(self)
            % Start the simulation run.
            %
            % Output:
            %   output - matrix of output values (as defined by model.out), prepended by time stamps
            
            % Check some simulation parameters
            self.checkSimParams();

            % The first set of output values at t = t_start:
            ut = self.input.u(self.t_start);
            initval = self.model.out(self.initial, ut, self.t_start);

            % To improve efficiency, preallocate the buffer for output:
            [OROWS, ~] = size(initval);
            output = self.preallocateOutput(OROWS+1);

            % Fill the 'output' with the initial value:
            output(:, 1) = [self.t_start; initval];

            s = self.initial;
            % Current index within 'output'
            idx = 2;
            for t = self.t_start : self.t_step : self.t_stop-self.t_step,
                sd = self.model.deriv(s, ut, t);
                s = s + sd * self.t_step;

                % Past this point, s represents states at the next point in time, i.e. at t+t_step.
                % This should be kept in mind when calculating output values and appliyng their time stamp.
                ut = self.input.u(t+self.t_step);
                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];

                % update 'idx'
                idx = idx+1;
            end % for
        end % function
    end % methods
    
end % classdef
