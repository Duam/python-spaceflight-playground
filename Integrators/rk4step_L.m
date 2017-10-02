function Lnext = rk4step_L (L, Lk, x, u, h)
% Implements one runge-kutta 4 step
% L: Stage to be integrated
% Lk: Current stage cost
% x: Current state
% u: Current control
% h: Stepsize

    k1 = L(x, u);
    k2 = L(x + h/2 * k1, u);
    k3 = L(x + h/2 * k2, u);
    k4 = L(x + h * k3, u);
    
    Lnext = Lk + h/6 * (k1 + 2*k2 + 2*k3 + k4);
