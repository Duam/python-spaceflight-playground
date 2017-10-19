clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D_cart.m

%% --------------------- Parameters --------------------------
% -- Simulation parameters --
T = 600;                % Simulation time in seconds
N = 100;                % Number of samples
DT = T/N;               % Step size

% -- Parameters of target orbit --
% Target altitude in km
h = 20;
R_orbit = 10^-3 * R + h; 
% Orbital angular velocity
angVel_T = sqrt(mu/(10^3 * R_orbit)^3); 

% -- Initial state --
x0 = [0; R; 0; 0; m0];

%% ------------------------ System model ----------------------------------
nx = 5;
nu = 2;
x = MX.sym('x', nx, 1);
u = MX.sym('u', nu, 1);
ode = Function('ode', {x,u}, {xdot(x,u)});
%% ---------------------- Stage cost --------------------------------------
L = u.' * u;
L = Function('L', {x,u}, {L});

%% ----------------------- Integrators ------------------------------------
nn = 10;   % Integrator steps per step
h = DT/nn; % Step size of integrator step

% Integration of the ODE
Xk = x;
for i=1:nn
    Xk = rk4step_ode(ode, Xk, u, h);
end
ode_d = Function('ode_d', {x,u}, {Xk}, {'x','u'}, {'xk'});

% Integration of the stage costs
Lk = 0;
for i=1:nn
    Lk = rk4step_L(L, Lk, x, u, h);
end
L_d = Function('L', {x,u}, {Lk}, {'x','u'}, {'Lk'});

%% ----------------------- Initial guess ----------------------------------
% -- Load precomputed trajectory as initial guess --
load('main/guess_PR2D_cart.mat');
x_g = sol.X;
u_g = sol.U;

%% Optimization variables
U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N);

%% Define OCP
% Start with empty NLP
w   = {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;
g   = {};
lbg = [];
ubg = [];

% Formulate NLP
Xk = x0;
for k=0:N-1
    % NLP variable for control
    Uk  = MX.sym(['U_' num2str(k)], nu, 1); 
    w   = {w{:}, Uk};
    lbw = [lbw; -inf; -inf];
    ubw = [ubw; inf;   inf];
    w0  = [w0; u_g(:,k+1)];

    % Circle constraints on the controls
    g = {g{:}, Uk.' * Uk};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    
    % Integrate till the end of the interval
    Xk_end = ode_d(Xk, Uk);
    J      = J + L_d(Xk, Uk);
    
    % New NLP variable for state
    Xk  = MX.sym(['X_' num2str(k+1)], nx);
    w   = {w{:}, Xk};
    lbw = [lbw; -inf; ...
                -inf; ...
                -inf; ...
                -inf; ...
                m_e ];
    ubw = [ubw; inf; ...
                inf; ...
                inf ; ...
                inf ; ...
                m0 ];
    w0  = [w0; x_g(:,k+1)];
    
    % Equality constraint to match intervals
    g   = {g{:}, Xk_end-Xk};
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];    
end

% Intermediate quantities
p = Xk_end(1:2);
v = Xk_end(3:4);
dist      = sqrt(p.' * p);
dist_err  = dist - R_orbit^2;
angVel     = (p(1)*v(2) - p(2)*v(1)) / dist;
angVel_err = angVel - angVel_T;
radVel     = (p(1)*v(1) + p(2)*v(2)) / dist;
radVel_err = radVel;

% End cost
J   = J - Xk_end(5).' * Xk_end(5);

% Terminal constraint on the altitude
g   = {g{:}, dist_err};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Terminal constraint on the angular velocity
g = {g{:}, angVel_err};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Terminal constraint on radial velocity
g = {g{:}, radVel_err};
lbg = [lbg; 0];
ubg = [ubg; 0];


% Clean up and reorder
g = vertcat(g{:});
w = vertcat(w{:});

% Create an NLP solver
nlp = struct('f', J, 'x', w, 'g', g);
options = struct;
options.ipopt.max_iter = 5000;
%options.ipopt.check_derivatives_for_naninf = 'yes';
options.ipopt.print_info_string = 'yes';
%options.ipopt.derivative_test = 'second-order';
solver = nlpsol('solver', 'ipopt', nlp, options);

% Solve the NLP
sol = solver('x0', w0, ...
             'lbx', lbw, ...
             'ubx', ubw, ...
             'lbg', lbg, ...
             'ubg', ubg);
w_opt = full(sol.x);
u_star = [w_opt(1:7:end), ...
          w_opt(2:7:end)].';
x_star = [w_opt(3:7:end), ...
          w_opt(4:7:end), ...
          w_opt(5:7:end), ...
          w_opt(6:7:end), ...
          w_opt(7:7:end)].';
x_star = [x0, x_star];
