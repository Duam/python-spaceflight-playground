clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D_cart.m
load('main/guess_PR2D_cart_tmp.mat');

%% --------------------- Parameters --------------------------
% -- Simulation parameters --
T = sol.param.T;                % Simulation time in seconds
N = sol.param.N;                % Number of samples
DT = T/N;               % Step size

% -- Parameters of target orbit --
% Target altitude in m
h_T = 20;
R_orbit = R + 10^3 * h_T; 
% Orbital angular velocity
angVel_T = 10^6 * sqrt(mu/R_orbit^3); 

% -- Initial state --
x0 = sol.x0;

%% ------------------------ System model ----------------------------------
nx = 5;
nu = 2;
x = MX.sym('x', nx, 1);
u = MX.sym('u', nu, 1);
ode = Function('ode', {x,u}, {xdot(x,u)});
%% ---------------------- Stage cost --------------------------------------
L = u.' * u; % Penalize controls
% L = L + (x(1:2).' * x(1:2) - 10^-3 * R_orbit); % Reach altitude
L = Function('L', {x,u}, {L});

%% ----------------------- Integrators ------------------------------------
nn = 10;   % Integrator steps per step
h_T = DT/nn; % Step size of integrator step

% Integration of the ODE
Xk = x;
for i=1:nn
    Xk = rk4step_ode(ode, Xk, u, h_T);
end
ode_d = Function('ode_d', {x,u}, {Xk}, {'x','u'}, {'xk'});

% Integration of the stage costs
Lk = 0;
for i=1:nn
    Lk = rk4step_L(L, Lk, x, u, h_T);
end
L_d = Function('L', {x,u}, {Lk}, {'x','u'}, {'Lk'});

%% ----------------------- Initial guess ----------------------------------
% -- Load precomputed trajectory as initial guess --
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
xend = scale.* Xk_end(:,end);
p = xend(1:2);
v = xend(3:4);
dist2 = p.' * p;
v_rad2 = (v.' * p)^2 / dist2;
v_tra2 = v.' * v - v_rad2;
omega2 = 10^12 * v_tra2 / dist2;

dist_err  = dist2 - R_orbit^2;
omega_err = angVel_T^2 - omega2;
v_rad_err = v_rad2;

% End cost
J   = J - Xk_end(5).' * Xk_end(5);

% Terminal constraint on the altitude
g   = {g{:}, dist_err};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Terminal constraint on the angular velocity
g = {g{:}, omega_err};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Terminal constraint on radial velocity
g = {g{:}, v_rad_err};
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
      
% Save the solution as a variable
param = struct('T', T, ...
               'N', N, ...
               'hT', h_T, ...
               'thetaDotT', angVel_T, ...
               'coordSys', 'cart');
sol = struct('x0', x0, ...
             'X', x_star, ...
             'U', u_star, ...
             'param', param);
             
