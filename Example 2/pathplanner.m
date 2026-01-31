function [path] = pathplanner(m1,m2,l,g,N2,inputConstraint)
% PATHPLANNER Generates an optimal trajectory for a cart-pendulum system
% using the Yop (CasADi) toolbox.
%
%   Inputs:
%   m1, m2, l, g    : Physical parameters
%   N2              : Number of control intervals for the refined solution
%   inputConstraint : [min, max] control input limits
%
%   DEPENDENCIES AND CITATIONS:
%   This function requires the Yop toolbox and CasADi.
%
%   1. Yop Toolbox
%      URL: https://github.com/yoptimalcontrol/yop
%      Ref: Leek, Viktor. "An optimal control toolbox for MATLAB based on CasADi." (2016).
%
%   2. CasADi
%      URL: https://web.casadi.org/
%      Ref: Andersson, Joel AE, Joris Gillis, Greg Horn, James B. Rawlings, and Moritz Diehl. 
%           "CasADi: a software framework for nonlinear optimization and optimal control." 
%           Mathematical Programming Computation 11, no. 1 (2019): 1-36.
%
%   3. IPOPT (Solver used by CasADi)
%      Ref: Wächter, Andreas, and Lorenz T. Biegler. "On the implementation of an 
%           interior-point filter line-search algorithm for large-scale nonlinear programming." 
%           Mathematical programming 106, no. 1 (2006): 25-57.

%% Cart-pendulum model
sys = YopSystem('states', 4, 'controls', 1);

% Symbolic variables
t = sys.t;
x = sys.x;
u = sys.u;
 
% Dynamics (Cart-Pendulum)
% x1: cart position, x2: cart velocity
% x3: pendulum angle, x4: pendulum angular velocity

xdot = [1;...
        (l*m2*sin(x(3))*(x(4))^2+u(1)+m2*g*cos(x(3))*sin(x(3)))/(m1+m2*(1-(cos(x(3)))^2));...
        x(4);...
        -(l*m2*cos(x(3))*sin(x(3))*(x(4))^2+u(1)*cos(x(3))+(m1+m2)*g*sin(x(3)))/(l*m1+l*m2*(1-(cos(x(3)))^2))];

sys.set('ode',xdot)

%% Optimal Control Problem (OCP) Formulation
ocp = YopOcp();

% Objective Function
ocp.min({t_f(x(1))})

% Constraints setup
ocp.st(... 
    'systems', sys, ...
    ... % Initial conditions
    {  0 '==' t_0(x(1))}, ... 
    {  0 '==' t_0(t)}, ...
    {  0 '==' t_0(x(2))}, ...
    {  0 '==' t_0(x(3))}, ...
    {  0 '==' t_0(x(4))}, ...
    ... % Terminal conditions
    {  3 '==' t_f(t)}, ...    
    {  0 '==' t_f(x(2))}, ...
    {  0 '==' t_f(x(3))}, ...
    {  0 '==' t_f(x(4))}, ...
    ... % Path constraints
    {  0 '<=' x(1) '<=' inf}, ... % Positive position only
    {  0 '<=' x(2) '<=' inf}, ... % Positive velocity only
    {  -0.15 '<=' x(3) '<=' 0.15}, ... % Angle limits
    { inputConstraint(1) '<=' u(1) '<=' inputConstraint(2)});

%% Solver executions

% Step 1: Solve on a coarse grid for warm start
N1 = 10;

 sol = ocp.solve( ...
     'controlIntervals', N1, ...
     'polynomialDegree', 1);

% Step 2: Solution using step 1 as warm start
 sol = ocp.solve( ...
     'controlIntervals', N2, ...
     'polynomialDegree', 1, ...
     'initialGuess',YopWarmstart(sol,sys,N2));
 
 %% Output Formatting
path = struct(...
    'Independent',sol.NumericalResults.Independent,...
    'State',sol.NumericalResults.State,...
    'Input',sol.NumericalResults.Control);

end

