%% ---------------- Parameters (All in kg, km, s) -------------------------
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

%% Rescaled ODE
% Original ODE: m, m/s, kg
% Scaled ODE: km, km/s, kg
scale = [10^3; 10^3; 10^3; 10^3; 1];
unscale = scale.^-1;

% Scale ODE
xdot = @(x,u) unscale .* ode(scale .* x, u, M);

%% ----------------- System ODE -------------------------------------------
% Coordinate frame origin is center of moon (Polar coordinates)
% nx = 5
% nu = 2
% x(1) -- X Position
% x(2) -- Y Position
% x(3) -- X Velocity
% x(4) -- Y Velocity
% x(5) -- Mass
% u(1) -- Input force in X
% u(2) -- Input force in Y

function [xdot] = ode(x, u, M)
    G  = 6.67408 * 10^-11;          % Grav. const. (m^3 * kg^-1 * s^-2)
    k_m   =      10^2;              % Fuel consumption coeff. (kg * s^-1)
    u_bar = 30 * 10^4;              % Max. thrust (N)

    
    p = x(1:2);
    m = x(5);
    
    coeff_grav = G*M / (p.' * p)^1.5;
    
    xdot = [x(3); ...
            x(4); ...
            u(1)*u_bar/m - x(1) * coeff_grav; ...
            u(2)*u_bar/m - x(2) * coeff_grav; ...
            -k_m * (u.' * u)];
end

