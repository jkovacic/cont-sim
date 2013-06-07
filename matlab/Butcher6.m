classdef Butcher6 < IntegRkAb
    % A differential equation solver, based on the 6th order Butcher's method.
    % The method is described at:
    % http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_butcher.c
    
    methods
        function obj = Butcher6(model, input, t_start, t_stop, t_step, init_condition)
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
            obj.A = [ 0, 0, 0, 0, 0, 0, 0; ... 
                      1/3, 0, 0, 0, 0, 0, 0; ... 
                      0, 2/3, 0, 0, 0, 0, 0; ...
                      1/12, 1/3, -1/12, 0, 0, 0, 0; ...
                      -0.0625,  1.125, -0.1875, -0.375, 0, 0, 0; ...
                      0, 1.125, -0.375, -0.75, 0.5, 0, 0; ...
                      9/44, -9/11, 63/44, 18/11, -16/11, 0, 0 ];
            obj.B = [11, 0, 81, 81, -32, -32, 11] / 120;
            obj.C = [0, 1/3, 2/3, 1/3, 0.5, 0.5, 1]';
            
        end % function
    
    end % methods
    
end % classdef
