classdef AB2 < IntegAdamsBashforthAb
    % A differential equation solver, based on the 
    % 1-step Adams - Bashforth - Moulton predictor - corrector method (actually the 
    % corrected Euler method). The method is described at:
    % http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
    % and
    % http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods

    
    methods
        function obj = AB2(model, input, t_start, t_stop, t_step, init_condition)
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
            
            obj = obj@IntegAdamsBashforthAb(model, input, t_start, t_stop, t_step, init_condition);
            
            % Coefficients for the 2-step Adams - Bashforth method:
            obj.ab_coef = [1.5, -0.5];
            
        end % function
    
    end % methods
    
end % classdef
