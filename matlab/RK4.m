classdef RK4 < IntegRkAb
    % A differential equation solver, based on the most common numerical method, 
    % i.e. the 4th order Runge - Kutta method. The method is described at:
    % http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Common_fourth-order_Runge.E2.80.93Kutta_method
    
    methods
        function obj = RK4(model, input, t_start, t_stop, t_step, init_condition)
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
            obj.A = [0, 0, 0, 0; 0.5, 0, 0, 0; 0, 0.5, 0, 0; 0, 0, 1, 0];
            obj.B = [1, 2, 2, 1] / 6;
            obj.C = [0, 0.5, 0.5, 1]';
            
        end % function
    
    end % methods
    
end % classdef
