#!/usr/bin/python3

#TODO: Is casadi needed here? I might need it when using
#      the rk4 step in an OCP
#import casadi as cas

##
# @brief One runge-kutta 4 step of an ordinary
#        differential equation xdot = f(x,u)
# @param f The time-continuous ODE to be integrated
# @param x The current state
# @param u The current control
# @param h The step size
# @return The state at the next timestep
##
def rk4step_ode(f, x, u, h):
    k1 = f(x,u)
    k2 = f(x + h/2 * k1, u)
    k3 = f(x + h/2 * k2, u)
    k4 = f(x + h * k3, u)
    
    return x + h/6 * (k1 + 2*k2 + 2*k3 + k4)

##
# @brief One runge-kutta 4 step of the lagrange term
#        in an optimal control problem
# @param L The time-continuous lagrange term
# @param Lk The current stage cost
# @param x The current state
# @param u The current control
# @param h The step size
# @return The next stage cost
##
def rk4step_L(L, Lk, x, u, h):
    k1 = L(x, u)
    k2 = L(x + h/2 * k1, u)
    k3 = L(x + h/2 * k2, u)
    k4 = L(x + h * k3, u)
    
    return Lk + h/6 * (k1 + 2*k2 + 2*k3 + k4)