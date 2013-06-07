classdef Gill4 < IntegRkAb
    % A differential equation solver, based on the Gill method.
    % The method is described at:
    % http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_gill.c
    
    
    methods
        function obj = Gill4(model, input, t_start, t_stop, t_step, init_condition)
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
            
            s2  = sqrt(2);

            % Elements of the Butcher tableau:
            obj.A = [0, 0, 0, 0; 1, 0, 0, 0; (-1+s2), (2-s2), 0, 0; 0, -s2, (2+s2), 0] / 2;
            obj.B = [1, (2-s2), (2+s2), 1] / 6;
            obj.C = [0, 0.5, 0.5, 1]';
            
        end % function
    
    end % methods
    
end % classdef
