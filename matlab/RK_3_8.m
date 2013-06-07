classdef RK_3_8 < IntegRkAb
    % A differential equation solver, based on the Runge - Kutta 3/8 method.
    % The method is described at:
    % http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    
    methods
        function obj = RK_3_8(model, input, t_start, t_stop, t_step, init_condition)
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
            obj.A = [0, 0, 0, 0; 1, 0, 0, 0; -1, 3, 0, 0; 3, -3, 3, 0] / 3;
            obj.B = [0.125, 0.375, 0.375, 0.125];
            obj.C = [0, 1, 2, 3]' / 3;
            
        end % function
    
    end % methods
    
end % classdef
