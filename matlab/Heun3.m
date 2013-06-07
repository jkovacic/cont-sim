classdef Heun3 < IntegRkAb
    % A differential equation solver, based on the 3rd order Heun method
    % (similar to the 3rd order Runge - Kutta method). The method is described at:
    % http://de.wikipedia.org/wiki/Runge-Kutta-Verfahren#Beispiele
    
    methods
        function obj = Heun3(model, input, t_start, t_stop, t_step, init_condition)
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
            obj.A = [0, 0, 0; 1, 0, 0; 0, 2, 0] / 3;
            obj.B = [1, 0, 3] / 4;
            obj.C = [0, 1, 2]' / 3;
            
        end % function
    
    end % methods
    
end % classdef
