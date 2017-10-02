function xnext = rk4step_ode (f, x, u, h)
% Implements one runge-kutta 4 step
% f: ODE to be integrated
% x: Current state
% u: Current control
% h: Stepsize

    k1 = f(x, u);
    k2 = f(x + h/2 * k1, u);
    k3 = f(x + h/2 * k2, u);
    k4 = f(x + h * k3, u);
    
    xnext = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
