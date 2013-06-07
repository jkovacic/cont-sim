classdef (Abstract) IntegAb
    % An abstract class that serves as a base for all derived differential equation
    % solving classes.
    %
    % It includes a declaration of equired properties and implementation of
    % some methods (mostly parameter checking), common to all derived
    % classes.
    
    properties (Access = protected)
        model;      % Model of a system (must be an implementation of IModel)
        input;      % External input(s) to the system (must be an implementation of IInput)
        
        t_start;    % Start time of the simulation run
        t_stop;     % Stop time of the simulation run
        t_step;     % Fixed time step
        
        initial;    % Vector of model's states' initial conditions
    end % properties
    
    methods (Access = protected)
        function obj = IntegAb(model, input, t_start, t_stop, t_step, init_condition)
            % A constructor that sets common simulation parameters.
            % It must be called from constructors of all derived classes.
            %
            % Input:
            %   model - state space model of the simulated system (must be an implementation of IModel)
            %   input - External input(s) to the system (must be an implementation of IInput), 
            %           if not required, InputNan can be passed
            %   t_start - start time of the simulation run
            %   t_stop - stop time of the simulation run
            %   t_step - fixed time step
            %   init_condition - initial condition, values of states at t=t_start
            
            obj.t_start = t_start;
            obj.t_stop = t_stop;
            obj.t_step = t_step;
            obj.initial = init_condition;
            
            obj.model = model;
            obj.input = input;
        end % function
        
        function checkSimParams(self)
            % Checks validity of simulation cycle parameters, e.g., all simulation run 
            % parameters must be scalars, t_stop must be greater than t_start,
            % t_step must be positive, model must be derivod from IModel, input must be 
            % derived from IInput, etc. If any parameter is invalid, a short
            % description of the problem will be displayed and an error will be thrown.
            
            if ( numel(self.t_start) ~= 1 )
                error ('Start time must be a scalar');
            end %if

            if ( numel(self.t_stop) ~= 1 )
                error ('Stop time must be a scalar');
            end %if

            if ( numel(self.t_step) ~= 1 )
                error ('Time step must be a scalar');
            end %if

            if ( self.t_stop < self.t_start )
                error ('Stop time must be greater than start time');
            end %if

            if ( self.t_step <= 0 )
                error('t_step must be positive');
            end % if
            
            if ( 0 == isa(self.model, 'IModel') )
                error('Model object must be derived from IModel');
            end %if
            
            if ( 0 == isa(self.input, 'IInput') )
                error('Model object must be derived from IInput');
            end %if
            
        end % function
        
        
        function output = preallocateOutput(self, rows)
            % Preallocates the buffer for output of fixed step algorithms
            % for improved efficiency.
            %
            % Input:
            %   rows - number of output's rows
            % Output:
            %   output - a [rows, floor((t_stop-t_start)/t_step)+1] matrix, filled with Nan 
            
            
            % It is assumed that regularity of parameters is checked 
            % by the caller beforehand

            output = NaN(rows, floor((self.t_stop-self.t_start)/self.t_step)+1);
        end %function
        
    end % methods
    
    methods
        % An abstract method (must be implemented by derived classes) that
        % runs a simulation cycle from t_start to t_stop and returns
        % model's outputs in 'output'
        %
        % Output:
        %   output - matrix of output values, prepended by time stamps
        output = run(self)
    end % methods
    
end % classdef
