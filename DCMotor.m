% Constants for Maxon DCX 22 S Motor and Blinds System
m = 2;                % kg, mass of blinds
R = 0.05;             % m, radius of blinds
J = 0.01;             % kg*m^2, moment of inertia of the blinds and motor system
Kt = 0.019;           % Nm/A, torque constant for DCX 22 S
Ke = 0.019;           % V/(rad/s), back EMF constant for DCX 22 S
R_m = 1.1;            % Ohms, motor resistance for DCX 22 S
L_m = 0.015;          % H, motor inductance for DCX 22 S
tau_friction = 0.05;  % Nm, friction torque in the system
omega_desired = 50;   % rad/s, desired motor speed

% PID Gains for controlling motor speed
Kp = 8;               % Proportional gain
Ki = 1;               % Integral gain
Kd = 1;               % Derivative gain

% Initial Conditions: [theta; omega; I; integral_error; prev_error]
init_conditions = [0; 0; 0.1; 0; 0];  % Initial motor state

% Simulation Time
tspan = [0 5]; % seconds

% Set ODE options for improved accuracy
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Run Simulation
[t, y] = ode45(@(t, y) blinds_dynamics(t, y, omega_desired, Kp, Ki, Kd, J, Kt, Ke, R_m, L_m, tau_friction), tspan, init_conditions, options);

% Calculate control voltage and speed error for each time step
control_voltage = zeros(size(t));    % Preallocate control voltage array
speed_error = zeros(size(t));        % Preallocate speed error array
for i = 1:length(t)
    % Extract omega (speed) from y
    omega = y(i, 2);
    
    % Calculate error and control voltage
    error = omega_desired - omega;
    speed_error(i) = error;
    integral_error = y(i, 4);
    derivative_error = (error - y(i, 5)) / 0.01; % approximate time step
    control_voltage(i) = Kp * error + Ki * integral_error + Kd * derivative_error;
end

% Plot Results
figure;
subplot(4,1,1);
plot(t, y(:,1));
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Blinds Position');

subplot(4,1,2);
plot(t, y(:,2));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Motor Speed');

subplot(4,1,3);
plot(t, control_voltage);
xlabel('Time (s)');
ylabel('Control Voltage (V)');
title('Control Voltage');

subplot(4,1,4);
plot(t, speed_error);
xlabel('Time (s)');
ylabel('Speed Error (rad/s)');
title('Speed Error');

% Dynamics Function
function dydt = blinds_dynamics(t, y, omega_desired, Kp, Ki, Kd, J, Kt, Ke, R_m, L_m, tau_friction)
    % Unpack Variables
    theta = y(1);       % Angular position
    omega = y(2);       % Angular velocity
    I = y(3);           % Motor current
    integral_error = y(4); % Integral of error for PID
    prev_error = y(5);  % Previous error for PID derivative
    
    % PID Control for Voltage V
    error = omega_desired - omega;
    integral_error = integral_error + error * 0.01; % Approximate integration step (0.01s)
    derivative_error = (error - prev_error) / 0.01; % Approximate derivative step (0.01s)
    
    V = Kp * error + Ki * integral_error + Kd * derivative_error;
    
    % Dynamics Equations
    % Electrical: L * dI/dt = V - R*I - Ke * omega
    dIdt = (V - R_m * I - Ke * omega) / L_m;
    
    % Mechanical: J * domega/dt = Kt * I - tau_friction
    domega_dt = (Kt * I - tau_friction) / J;
    
    % dTheta/dt = omega
    dtheta_dt = omega;
    
    % Resulting Differential Equations
    dydt = [dtheta_dt; domega_dt; dIdt; error; error];
end
