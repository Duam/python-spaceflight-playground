%% ---------------- Parameters (All in kg, km, s) --------------------------
% -- Physical parameters --
G  = 6.67408 * 10^-11;          % Grav. const. (m^3 * kg^-1 * s^-2)
M  = 7.348 * 10^22;             % Moon mass (kg)
mu = G * M;                     % Std. grav. param. (km^3 * s^-2)
R  = 1737.5 * 10^3;             % Moon radius (m)

% -- Rocket parameters --
m0    = 20 * 10^3;              % Init. rocket mass (kg)
m_e   =  1 * 10^3;              % Empty rocket mass (kg)
u_bar = 30 * 10^4;              % Max. thrust (N)
k_m   =      10^2;              % Fuel consumption coeff. (kg * s^-1)

%% ----------------- System ODE -------------------------------------------
% Coordinate frame origin is center of moon (Polar coordinates)
% nx = 5
% nu = 2
% x(1) -- Distance from origin
% x(2) -- Angle
% x(3) -- Velocity
% x(4) -- Angular velocity
% x(5) -- Mass
% xdot = @(x,u) [
%     x(3); ...
%     x(4); ...
%     u(1)*u_bar/x(5)  - mu/x(1)^2 + x(1)*x(4)^2; ...
%     u(2)*u_bar/(x(5)*x(1)) - 2*x(3)*x(4)/x(1) ; ...
%     - k_m * (u(1)^2 + u(2)^2)
%     ];

% System with state 1 = altitude
xdot = @(x,u) [
    x(3); ...
    x(4); ...
    u(1)*u_bar/x(5)  - mu/(x(1)+R)^2 + (x(1)+R)*x(4)^2; ...
    u(2)*u_bar/(x(5)*(x(1)+R)) - 2*x(3)*x(4)/(x(1)+R) ; ...
    - k_m * (u(1)^2 + u(2)^2)
    ];

% Scaling factor
% km becomes m
% nrad becomes rad
% kg stays kg
scale = [10^3; 10^-9; 10^3; 10^-9; 1];
unscale = scale.^-1;

% Scale ODE
xdot = @(x,u) unscale .* xdot(scale .* x, u);
