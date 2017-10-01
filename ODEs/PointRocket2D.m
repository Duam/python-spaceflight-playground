%% ---------------- Parameters (All in kg, m, s) --------------------------
% -- Physical parameters --
G = 6.67408 * 10^-11;   % Gravitational constant
M = 7.348 * 10^22;      % Moon mass
mu = G * M;             % Standard gravitational parameter
R = 1.7375 * 10^6;      % Moon radius

% -- Rocket parameters --
m0   = 20 * 10^3;      % Initial mass of rocket
m_e   =  1 * 10^3;      % Fuel-empty mass of rocket
u_bar = 30 * 10^4;      % Maximum radial thrust
k_m   =      10^2;      % Fuel consumption coefficient

%% ----------------- System ODE -------------------------------------------
% Coordinate frame origin is center of moon (Polar coordinates)
% nx = 5
% nu = 2
% x(1) -- Distance from origin
% x(2) -- Angle
% x(3) -- Velocity
% x(4) -- Angular velocity
% x(5) -- Mass
xdot = @(x,u) [
    x(3); ...
    x(4); ...
    u(1)*u_bar/x(5)  - mu/x(1)^2 + x(1)*x(4)^2; ...
    u(2)*u_bar/(x(5)*x(1)) - 2*x(3)*x(4)/x(1) ; ...
    - k_m * (u(1)^2 + u(2)^2)
    ];