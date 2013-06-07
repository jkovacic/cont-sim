classdef RK2 < IntegRkAb
    % A differential equation solver, based on the the 2nd order Runge - Kutta method
    % (alpha=1/2. a.k.a "the midpoint method").
    % The method is described at:
    % http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Second-order_methods_with_two_stages
    
    methods
        function obj = RK2(model, input, t_start, t_stop, t_step, init_condition)
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
            [obj.A, obj.B, obj.C] = IntegRkAb.rk2Butcher(0.5);
        end % function
    
    end % methods
    
end % classdef
