classdef Nystrom5 < IntegRkAb
    % A differential equation solver, based on the 5th order Nystroem's method.
    % The method is described at:
    % http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_nystrom.c
    
    methods
        function obj = Nystrom5(model, input, t_start, t_stop, t_step, init_condition)
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
            obj.A = [ 0, 0, 0, 0, 0, 0; ...
                      1/3, 0, 0, 0, 0, 0; ...
                      0.16, 0.24, 0, 0, 0, 0; ...
                      0.25, -3, 15/4, 0, 0, 0; ...
                      2/27, 10/9, -50/81, 8/81, 0, 0; ...
                      0.08, 0.48, 2/15, 8/75, 0, 0];
            obj.B = [23, 0, 125, 0, -81, 125] / 192;
            obj.C = [0, 1/3, 0.4, 1, 2/3, 0.8]';
            
        end % function
    
    end % methods
    
end % classdef
