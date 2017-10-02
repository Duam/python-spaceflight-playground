clear all
import casadi.*

%% Load ODE and parameters
% Current working directory must be "ControlledRocket"
run ODEs/PointRocket2D.m

%% ---------------- Parameters (All in kg, m, s) --------------------------
% -- Simulation parameters --
T = 600;                % Simulation time in seconds
N = 100;                % Number of samples
DT = T/N;               % Step size

% -- Parameters of target orbit --
R_T = R + 20000;        % Target radius

% -- Initial state --
x0 = [R; pi/2; 0; 0; m0];   % Rocket starts from "top" of the moon

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
winit = load('main/guess_PR2D.mat');
x_g = winit.winit_PR2D.State;
u_g = winit.winit_PR2D.Control;

% Next: scale the ode

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
    lbw = [lbw; R; ...
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
dR_T        = Xk_end(1) - R_T;
dThetaDot_T = Xk_end(4) - sqrt(mu/R_T^3);

% End cost
J   = J - Xk_end(5).' * Xk_end(5);

% Terminal constraint on the angular velocity
g   = {g{:}, dThetaDot_T};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Terminal constraint on the altitude
g = {g{:}, dR_T};
lbg = [lbg; 0];
ubg = [ubg; 0];

% Clean up and reorder
g = vertcat(g{:});
w = vertcat(w{:});

% Create an NLP solver
nlp = struct('f', J, 'x', w, 'g', g);
options = struct;
options.ipopt.max_iter = 5000;
options.ipopt.check_derivatives_for_naninf = 'yes';
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
