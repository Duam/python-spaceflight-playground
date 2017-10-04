%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D.m

% Load results
% results = load('results/PR2D_sol1.mat');
% X = results.PR2D_sol1.X;
% U = results.PR2D_sol1.U;
results = load('results/PR2D_sol2.mat');
X = results.sol.X;
U = results.sol.U;


% Simulation parameters
T = 600;
N = 100;
DT = T/N;

% Target altitude
h_T = 20;


%% -- Plot --
tAxis = 0:DT:T;
% Radius of planet
r_0 = zeros(1,length(tAxis));
r_T = h_T * ones(1,length(tAxis));
% Angle
theta_0 = pi/2 * ones(1,length(tAxis));
% Angular velocity limits
thetaDot_0 = zeros(1,length(tAxis));
thetaDot_T = 10^6 * sqrt(G*M/(R+h_T)^3) * ones(1,length(tAxis));
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
grid on
% Plot angle
subplot(3,2,2);
hold on
plot(tAxis, theta_0, '--r');
plot(tAxis, X(2,:));
ylabel('$\theta$', 'interpreter', 'latex');
grid on
% Plot radial velocity
subplot(3,2,3);
hold on
plot(tAxis, zeros(1,N+1), '--r');
plot(tAxis, X(3,:));
ylabel('$\dot{r}$', 'interpreter', 'latex');
grid on
% Plot angular velocity
subplot(3,2,4);
hold on
plot(tAxis, thetaDot_0, '--r');
plot(tAxis, thetaDot_T, '--r');
plot(tAxis, X(4,:));
ylabel('$\dot{\theta}$', 'interpreter', 'latex');
grid on
% Plot mass
subplot(3,2,[5 6]);
hold on
plot(tAxis, m_0, '--r');
%plot(tAxis, m_T, '--r');
plot(tAxis, X(5,:));
ylabel('$m$', 'interpreter', 'latex');
grid on