clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D_cart.m

%% Parameters (All in kg, m, s)
% Simulation parameters
T = 600;
N = 100;  % Number of samples
DT = T/N; % Step size

% Initial state (km, km/s, kg)
x0 = unscale .* [0; R; 0; 0; m0];

%% System model
nx = 5;
nu = 2;
x = MX.sym('x', nx, 1);
u = MX.sym('u', nu, 1);

ode = Function('ode', {x,u}, {xdot(x,u)});

%% RK4-step for ode 
nn = 10;   % Integrator steps per step
h = DT/nn; % Step size of integrator step
Xk = x;
for i=1:nn % Oversampling
    Xk = rk4step_ode(ode, Xk, u, h);
end
% Discretized ode
ode_d = Function('ode_d', {x,u}, {Xk}, {'x','u'}, {'xk'});

%% Simulate system
% Controls
% nx_stop = 80 ;
% Ux = [zeros(1,(N-nx_stop)/2), -0.5 * ones(1,nx_stop), zeros(1,(N-nx_stop)/2)];
% ny_stop = 30;
% Uy = [0.125 * ones(1,ny_stop), zeros(1,N-ny_stop)];
% U = [Ux; Uy];
load('Results/PR2D_sol2_convertedToCart.mat', 'sol');
U = sol.U;

% Simulate
X(:,1) = x0;
for i=2:N+1
    X(:,i) = full(ode_d(X(:,i-1), U(:,i-1)));
end

% Altitude of target orbit
h_T = 20;
% Orbital angular velocity in microradians
angVel_T = 10^6 * sqrt(mu/(R+10^3*h_T)^3); 

% Radial velocity @ end (magnitude)
% v_rad = (v2 - v1) dot (r2 - r1)/|r2 - r1|
% v1 = 0, r1 = 0 (inertial frame)
x_end = scale.* X(:,end);
p = x_end(1:2);
v = x_end(3:4);
dist = sqrt(p.' * p); % must be == R + 10^3 * h
v_rad = (v.' * p) / dist; % must be == 0
% Transverse velocity @ end (magnitude)
% v_tra = sqrt(|v2 - v1|^2 - v_rad^2)
v_tra = sqrt(v.' * v - v_rad^2);
% Angular velocity @ end
omega = 10^6 * v_tra / dist; % must be == angVel_T

disp('Distance: ' + string(dist) + ', target: ' + string(R+10^3*h_T) + ', Diff: ' + string(R+10^3*h_T-dist));
disp('omega:    ' + string(omega) + ', target: ' + string(angVel_T) + ', Diff: ' + string(angVel_T - omega));
disp('RadVel:   ' + string(v_rad) + ', target: ' + string(0) + ', Diff: ' + string(v_rad));

% Save the solution as a variable
param = struct('T', T, ...
               'N', N, ...
               'hT', h_T, ...
               'thetaDotT', angVel_T, ...
               'coordSys', 'cart');
sol = struct('x0', x0, ...
             'X', X(:,2:end), ...
             'U', U, ...
             'param', param);
       
% Plot
PR2D_plotResults(sol);

save('main/guess_PR2D_cart_tmp.mat', 'sol');