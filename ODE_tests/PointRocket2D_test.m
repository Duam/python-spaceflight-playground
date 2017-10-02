clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D.m

%% Parameters (All in kg, m, s)
% Target orbit
T_orbit = 2*pi*sqrt(R^3 / mu) * 2;
R_T = R + 20000;

% Simulation parameters
T = 600;
N = 100;  % Number of samples
DT = T/N; % Step size

% Initial state
x0 = [R; pi/2; 0; 0; m0];

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
    Xk = rk4step(ode, Xk, u, h);
end
% Discretized ode
ode_d = Function('ode_d', {x,u}, {Xk}, {'x','u'}, {'xk'});

%% Simulate system
n_r_stop = 60 ;
U_rad = [0.1075 *ones(1,n_r_stop), zeros(1,N-n_r_stop)];
n_theta_stop = 85;
U_tan = [0.2 * ones(1, n_theta_stop), zeros(1,N-n_theta_stop)];
% U_tan = zeros(1,N);
U = [U_rad; U_tan];

X(:,1) = x0;
for i=2:N
    X(:,i) = full(ode_d(X(:,i-1), U(:,i-1)));
end

%% -- Plot --
tAxis = 0:DT:T-DT;
% Radius of planet
r_0 = R * ones(1,length(tAxis));
r_T = R_T * ones(1,length(tAxis));
% Angle
theta_0 = pi/2 * ones(1,length(tAxis));
% Angular velocity limits
thetaDot_0 = zeros(1,length(tAxis));
thetaDot_T = sqrt(G*M/R_T^3) * ones(1,length(tAxis));
% Mass limits
m_0 = m0 * ones(1,length(tAxis));
m_T = m_e * ones(1,length(tAxis));

figure(1);
clf
% Plot radius
subplot(3,2,1);
hold on
plot(tAxis, r_0, '--r');
plot(tAxis, r_T, '--r');
plot(tAxis, X(1,:));
ylabel('$r$', 'interpreter', 'latex');
% Plot angle
subplot(3,2,2);
hold on
plot(tAxis, theta_0, '--r');
plot(tAxis, X(2,:));
ylabel('$\theta$', 'interpreter', 'latex');
subplot(2,1,2);
% Plot radial velocity
subplot(3,2,3);
hold on
plot(tAxis, zeros(1,N), '--r');
plot(tAxis, X(3,:));
ylabel('$\dot{r}$', 'interpreter', 'latex');
% Plot angular velocity
subplot(3,2,4);
hold on
plot(tAxis, thetaDot_0, '--r');
plot(tAxis, thetaDot_T, '--r');
plot(tAxis, X(4,:));
ylabel('$\dot{\theta}$', 'interpreter', 'latex');
% Plot mass
subplot(3,2,[5 6]);
hold on
plot(tAxis, m_0, '--r');
plot(tAxis, m_T, '--r');
plot(tAxis, X(5,:));
ylabel('$m$', 'interpreter', 'latex');


