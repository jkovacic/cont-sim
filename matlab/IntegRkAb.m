classdef (Abstract) IntegRkAb < IntegAb
    % An abstract class with common functionality for all explicit Runge -
    % Kutta based methods for solving differential equations.
    %
    % All classes based on Runge - Kutta like algorithms should be derived
    % from this one.
    
    properties (Access=protected)
        A;  % matrix A of the Butcher tableau (its diagonal and upper triangle will never be read)
        B;  % matrix B of the Butcher tableau
        C;  % matrix C of the Butcher tableau (the first element will automatically be set to 0)
    end % properties
    
    methods (Access=protected)
        function obj = IntegRkAb(model, input, t_start, t_stop, t_step, init_condition)
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
            
            obj = obj@IntegAb(model, input, t_start, t_stop, t_step, init_condition);
        end % function
    end % methods
    
    methods (Access=private)
        
        function checkButcherElements(self)
            % Checks whether elements of the Butcher tableau (matrices A, B
            % and C) are valid. If any matrix has invalid dimensions or 
            % matrices' elements are inconsistent, an error message will be 
            % thrown. The first element of C and the
            % upper diagonal elements of A will be set to 0.
            %
            % For more details about elements of a Butcher tableau, see:
            % http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
            
            % Check validity of Butcher tableau elements:
            % A must be a square matrix:
            [N, columns] = size(self.A);
            if ( N ~= columns )
                error('A must be square matrix');
            end %if

            % B should be a 1xN vector, however any matrix with N elemnts will do
            if ( numel(self.B) ~= N )
                error('Invalid size of B');
            end %if

            % C should be a Nx1 vector, however any matrix with N elemnts will do
            if ( numel(self.C) ~= N )
                error('Invalid size of C');
            end %if

            % Some elements of C and A should be 0.
            % Instead of checking their values and reporting an error,
            % those elements will be set to zeros:
            self.C(1) = 0;
            for i = 1:N
                self.A(i, i:N) = zeros(1, N-i+1);
            end %for

            % Consistency of the Butcher tableau:
            % Note, due to limited accuracy of floating point arithmetics,
            % a bit larger threshold for "equality" is necessary for some methods, e.g. Ralston's or Verner's.
            crit = 7 * eps;
            for i = 1:N
                % This might be useful to determine the threshold
                % diff = sum(A(i, :)) - C(i);
                % printf("diff: %f    crit: %f\n", diff/eps, crit/eps);
                if ( abs( sum(self.A(i, :)) - self.C(i) ) > crit )
                    error('Inconsistent Butcher tableau');
                end %if
            end %for
        end % function
        
    end % methods
    
    
    methods (Static)
        function [A, B, C] = rk2Butcher(alpha)
            % Elements of a general Butcher tableau for 2nd order Runge - Kutta methods.
            % For more details, see:
            % http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Second-order_methods_with_two_stages
            %
            % As the method might be useful elsewhere outside of this class, 
            % it is implemented as a public static method
            %
            % Input:
            %   alpha - parameter of a 2nd order Runge - Kutta method
            %
            % Output:
            %   A - matrix A of the Butcher tableau
            %   B - vector B of the Butcher tableau
            %   C - vector C of the Butcher tableau
            
            if ( abs(alpha) <= eps )
                error('Alpha cannot be zero!');
            end %if

            A = zeros(2, 2);
            B = zeros(1, 2);
            C = zeros(2, 1);

            % for a slight improvement in speed, calculate 1/(2*alpha) only once:
            invalpha = 1/(2*alpha);

            A(2, 1) = alpha;
            C(2) = alpha;

            B(1) = 1 - invalpha;
            B(2) = invalpha;

            end %function
        
    end % methods
    
    methods
        function output = run(self)
            % Start the simulation run.
            % The common "engine" is used for all explicit Runge - Kutta
            % based methods.
            %
            % Output:
            %   output - matrix of output values (as defined by model.out), prepended by time stamps

            self.checkButcherElements();
            
            [N, ~] = size(self.A);
            
            % Check some simulation parameters
            self.checkSimParams();

            % Dimensions of a matrix of states. Even though it is recommended to be a vertical
            % vector (n x 1), the implementation is robust enough to handle other dimensions as well.
            [STATE_ROWS, STATE_COLS] = size(self.initial);

            % Check validity of the vector/matrix of states and coefficients:
            if ( STATE_COLS <=0 )
                error('Invalid dimensions of initial_condition.');
            end %if

            % Matrix to store values of kn
            K = zeros(STATE_ROWS, N*STATE_COLS);

            % The first set of output values at t = t_start:
            ut = self.input.u(self.t_start);

            initval = self.model.out(self.initial, ut, self.t_start);

            % To improve efficiency, preallocate the buffer for output:
            [NOUT, ~] = size(initval);
            output = self.preallocateOutput(1+NOUT);

            % Fill the 'output' with the initial value:
            output(:, 1) = [self.t_start; initval];


            % Start of the Runge - Kutta algorithm:

            s = self.initial;
            % Current index within 'output'
            idx = 2;
            for t = self.t_start : self.t_step : self.t_stop-self.t_step
                % The current value of s is used several times inside the for i loop.
                % However, to eliminate the need for another foor loop (after for i),
                % s can be updated inside the same for loop. As this would break the algorithm
                % for calculation of kn, stemp was introduced.
                stemp = s;
                for i = 1:N
                    rktemp = stemp;
                    for j = 1: (i-1)
                        rktemp = rktemp + self.A(i, j) * K(:, ((j-1)*STATE_COLS+1) : (j*STATE_COLS) );
                    end %for

                    % ki and its position in K:
                    ui = self.input.u(t+self.C(i)*self.t_step);
                    K(:, ((i-1)*STATE_COLS+1) : (i*STATE_COLS) ) = ...
                        self.t_step * self.model.deriv(rktemp, ui, t+self.C(i)*self.t_step);

                    s = s + self.B(i) * K(:, ((i-1)*STATE_COLS+1) : (i*STATE_COLS) );
                end %for

                % Past this point, s represents states at the next point in time, i.e. at t+t_step.
                % This should be kept in mind when calculating output values and applying their time stamp.
                ut = self.input.u(t+self.t_step);
                val = self.model.out(s, ut, t+self.t_step);
                output(:, idx) = [t+self.t_step; val];
                % update 'idx'
                idx = idx + 1;
            end % for

        end %function
    end % methods
    
end % classdef
