function w0 = YopWarmstart(sol,sys,N,x_k,u_sol,xp0)

t = sol.signal(sys.t);
x = sol.signal(sys.x);
u = sol.signal(sys.u);
u_resampled = interp1(t',u',0:1/N:1)'; % has to do with small scale
x_resampled = interp1(t',x',0:1/N:1)';
t_resampled = interp1(t',t',0:1/N:1);
W  = [t_resampled' x_resampled' u_resampled'];
% Create up-sampled as initial guess
w0 = YopInitialGuess( ...
    'signals', [sys.t; sys.x; sys.u], ...
    'signalValues',W');