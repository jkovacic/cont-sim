classdef Ralston4 < IntegRkAb
    % A differential equation solver, based on the 4th order Ralston's method.
    % The method is described at:
    % http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_ralston_4.c
    
    methods
        function obj = Ralston4(model, input, t_start, t_stop, t_step, init_condition)
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
            
            s5 = sqrt(5);

            % Elements of the Butcher tableau:
            obj.A = [0, 0, 0, 0; 0.4, 0, 0, 0; (-2889+1428*s5)/1024, (3785-1620*s5)/1024, 0, 0; ...
                     (-3365+2094*s5)/6040, (-975-3046*s5)/2552, (467040+203968*s5)/240845, 0];
            obj.B = [(263+24*s5)/1812, (125-1000*s5)/3828, 1024*(3346+1623*s5)/5924787, (30-4*s5)/123];
            obj.C = [0, 0.4, 0.875-0.1875*s5, 1]';
            
        end % function
    
    end % methods
    
end % classdef
