classdef MovablePendulum < IModel
    % A model of a pendulum on a movable support as derived at:
    % http://en.wikipedia.org/wiki/Lagrangian_mechanics#Pendulum_on_a_movable_support
    
    % The model's states have the following meanings:
    % s(1) = x [m]          - horizontal position of the pendulum support 
    % s(2) = xd [m/s]       - dx/dt, horizontal speed of the pendulum support
    % s(3) = theta [rad]    - inclination of the pendulum (relative to the vertical line)
    % s(4) = thetad [rad/s] - dtheta/dt, rotational speed of the pendulum
    
    properties (Constant)
        N_STATES = 4;   % Number of state variables
    end % properties

    properties
        % Model's parameters (set by the constructor). To facilitate
        % possible experiments with parametrization etc., all parameters
        % are public properties.
        
        M;  % [kg]    - mass of the movable body
        m;  % [kg]    - mass of the pendulum point
        l;  % [m]     - length of the pendulum
        g;  % [m/s^2] - gravitational acceleration
    end % properties
    
    methods
        
         function obj = MovablePendulum(param)
             % A constructor that sets values of parameters.
             %
             % Meanings of each parameter:
             %   param(1) = M [kg]    - mass of the movable body
             %   param(2) = m [kg]    - mass of the pendulum point
             %   param(3) = l [m]     - length of the pendulum
             %   param(4) = g [m/s^2] - gravitational acceleration
             %
             % Input:
             %   params - a vector of parameters (at least 4 numeric values)
     
             % Check that at least 4 parameters were provided
             % (note that this number has no relation with N_STATES):
             if ( 4 ~= length(param) )
                 error ('At least 4 elements expected in ''params''');
             end
             
             obj.M = param(1);
             obj.m = param(2);
             obj.l = param(3);
             obj.g = param(4);
         end % function
        
         function sd = deriv(self, s, ~, ~)
            % Calculates derivatives of model's state variables with
            % respect of time.
            %
            % Input:
            %   s - vector of values of states at the specified moment in time
            %   u - external input at the sepecified moment in time (not used by this method)
            %   t - moment in time (not used by this method)
            % Output:
            %   sd - vector of derivatives of states

            if ( MovablePendulum.N_STATES ~= numel(s) )
                error('Invalid number of states, expected %d', MovablePendulum.N_STATES);
            end %if

            % As 'x' itself is not used at the model, it is commented out:
            %x = s(1);        % [m]     - horizontal position of the pendulum support 
            xd = s(2);       % [m/s]   - dx/dt, horizontal speed of the pendulum support
            theta = s(3);    % [rad]   - inclination of the pendulum (relative to the vertical line)
            thetad = s(4);   % [rad/s] - dtheta/dt, rotational speed of the pendulum
            
            % The derived model is a set of two ordinary differential equations:
            %
            % (M+m)*xdd + m*l*cos(theta)*thetadd - m*l*sin(theta)*thetad^2 = 0
            % thetadd + xdd*cos(theta)/l + g*sin(theta)/l = 0                      (1)
            %
            % The highest order derivatives of both states (xdd and thetadd) must be expressed
            % as explicite functions of the lower level derivatives, inputs, time, params, etc.
            %
            % From the first equation, xdd can be derived as:
            % xdd = m*l*(thetad^2*sin(theta)-thetadd*cos(theta))/(m+M)             (2)
            %
            % and put into the the second equation which turns into:
            % thetadd + m*(thetad^2*sin(theta)*cos(theta)-thetadd*cos^2(theta))/(m+M) +g*sin(theta)/l = 0       (3)
            %
            % thetadd can be derived from it as:
            % thetadd = (-sin(theta)*(m*thetad^2*cos(theta)/(m+M))+g/l)/(1-m*cos^2(theta)/(m+M))                (4)
            %
            % Put the thetadd into (2) and xdd can be expressed as:
            % xdd = (m*l/(m+M))*(thetad^2*sin(theta)+(sin(theta)*cos(theta)*(m*thetad^2*cos(theta)/(m+M))+g/l)/ \\
            %   (1-m*cos^2(theta)/(m+M)))                                                                       (5)
            %
            % The equations (4) and (5) can be used to calculate sd.
            % Note that (1-m*cos^2(theta)/(m+M)) can never limit to zero unless M is negligible comparing to m.
            % As this is unlikely in typical situations, the division by zero is very unlikely
            % and will not be handled separately


            % Sine and cosine of theta are used often in the model. Since trigonometrical functions are 
            % computationally complex, store the theta's sine and cosine into temporary values:
            sth = sin(theta);
            cth = cos(theta);

            % m/(m+M) and g/l also appear often, so the values will be stored into separate variables:
            mrel = self.m/(self.m+self.M);
            gl = self.g/self.l;

            % Finally calculate numerical values of both second order derivatives as derived above:
            xdd = mrel*self.l*(thetad*thetad*sth + (sth*cth*(mrel*thetad*thetad*cth+gl)) / (1-mrel*cth*cth));
            thetadd = -sth*(mrel*thetad*thetad*cth+gl)/(1-mrel*cth*cth);

            % And appropriately fill the output vector of derivatives:
            sd = s;    % to preserve the dimensions and avoid errors due to addition of incompatible vectors/matrices
            sd(1) = xd;
            sd(2) = xdd;
            sd(3) = thetad;
            sd(4) = thetadd;

        end % function

        function output = out(~, s, ~, ~)
            % Calculates desired outputs from the vector of states.
            %
            % Input:
            %   s - vector of values of states at the specified moment in time
            %   u - external input at the sepecified moment in time (not used by this method)
            %   t - moment in time (not used by this method)
            % Output:
            %   out - vector of desired output values at the specified time
            %
            % The desired output vector consists of the following values:
            %   x [cm]      - horizontal position of the pendulum support 
            %   theta [deg] - inclination of the pendulum (relative to the vertical line)
            
            if ( MovablePendulum.N_STATES ~= numel(s) )
                error('Invalid number of states, expected %d', MovablePendulum.N_STATES);
            end %if

            % Note that x must be coverted from meters to centimeters (multiplying by 100)
            % and theta must be coverted from radians to degrees (multiplying by 180/pi).
            output = NaN(2,1);
            output(1) = s(1) * 100;
            output(2) = s(3) * 180/pi;
            
            % Could also be impleneted as:
            % out = [s(1)*100, s(3)*180/pi]';

        end % function
        
    end % methods
    
end % classdef
