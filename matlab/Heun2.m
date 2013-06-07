classdef Heun2 < IntegRkAb
    % A differential equation solver, based on the the 2nd order Heun method 
    % (similar to the 2nd Runge - Kutta method, alpha=1).
    % The method is described at:
    % http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Second-order_methods_with_two_stages
    % and
    % http://en.wikipedia.org/wiki/Heun%27s_method
    
    methods
        function obj = Heun2(model, input, t_start, t_stop, t_step, init_condition)
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
            
            obj = obj@IntegRkAb(model, input, t_start, t_stop, t_step, init_condition);
            
            % Elements of the Butcher tableau:
            [obj.A, obj.B, obj.C] = IntegRkAb.rk2Butcher(1);
        end % function
    
    end % methods
    
end % classdef
