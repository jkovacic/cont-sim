classdef ABM3 < IntegAdamsBashforthMoultonAb
    % A differential equation solver, based on the  
    % 3-step Adams - Bashforth - Moulton predictor - corrector method. The method is described at:
    % http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Bashforth_methods
    % and
    % http://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods
    
    methods
        function obj = ABM3(model, input, t_start, t_stop, t_step, init_condition)
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
            
            obj = obj@IntegAdamsBashforthMoultonAb(model, input, t_start, t_stop, t_step, init_condition);
            
            % Coefficients for the 3-step Adams - Bashforth method:
            obj.ab_coef = [23, -16, 5] / 12;
            
            % Coefficients for the 4-step Adams - Moulton method:
            obj.am_coef = [9, 19, -5, 1] / 24;
            
        end % function
    
    end % methods
    
end % classdef
