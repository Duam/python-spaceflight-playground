clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D.m

%% Parameters (All in kg, m, s)
% Altitude of target orbit
h_T = 20;
% Orbital angular velocity in microradians
angVel_T = 10^6 * sqrt(mu/(R+10^3*h_T)^3); 

% Simulation parameters
T = 600;
N = 100;  % Number of samples
DT = T/N; % Step size

% Initial state
x0 = [0; 0; 0; 0; m0];

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
n_r_stop = 60 ;
U_rad = [0.1075 *ones(1,n_r_stop), zeros(1,N-n_r_stop)];
n_theta_stop = 85;
U_tan = [0.2 * ones(1, n_theta_stop), zeros(1,N-n_theta_stop)];
U = [U_rad; U_tan];

% Simulate
X(:,1) = x0;
for i=2:N
    X(:,i) = full(ode_d(X(:,i-1), U(:,i-1)));
end

%% -- Plot --
tAxis = 0:DT:T-DT;
r_0        = x0(1) * ones(1,length(tAxis));
theta_0    = x0(2) * ones(1,length(tAxis));
rDot_0     = x0(3) * ones(1,length(tAxis));
thetaDot_0 = x0(4) * ones(1,length(tAxis));
m_0        = x0(5) * ones(1,length(tAxis));

r_T        = h_T * ones(1,length(tAxis));
thetaDot_T = angVel_T * ones(1,length(tAxis));
m_T        = m_e * ones(1,length(tAxis));

figure(1);
clf
% Plot radius
subplot(3,2,1);
hold on
plot(tAxis, r_0, '--r');
plot(tAxis, r_T, '--r');
plot(tAxis, X(1,:));
ylabel('$r [km]$', 'interpreter', 'latex');
% Plot angle
subplot(3,2,2);
hold on
plot(tAxis, theta_0, '--r');
plot(tAxis, X(2,:));
ylabel('$\theta [\mu rad]$', 'interpreter', 'latex');
subplot(2,1,2);
% Plot radial velocity
subplot(3,2,3);
hold on
plot(tAxis, rDot_0, '--r');
plot(tAxis, X(3,:));
ylabel('$\dot{r} [\frac{km}{s}]$', 'interpreter', 'latex');
% Plot angular velocity
subplot(3,2,4);
hold on
plot(tAxis, thetaDot_0, '--r');
plot(tAxis, thetaDot_T, '--r');
plot(tAxis, X(4,:));
ylabel('$\dot{\theta} [\frac{\mu rad}{s}]$', 'interpreter', 'latex');
% Plot mass
subplot(3,2,[5 6]);
hold on
plot(tAxis, m_0, '--r');
plot(tAxis, m_T, '--r');
plot(tAxis, X(5,:));
ylabel('$m [kg]$', 'interpreter', 'latex');


