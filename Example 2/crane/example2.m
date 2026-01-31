clear; close all; clc;
%% Obtaining the reference path
% System Physical Parameters
m1 = 1.2;       % Mass of cart
m2 = 0.6;       % Mass of pendulum
l  = 1.5;       % Length of pendulum
g  = 9.8067;    % Gravity

% Controller / Planner Settings
ctrl_intervals = 100;       % Number of control intervals (resolution)
ctrl_const     = [-2, 2];   % Input constraints [min, max]
createPath     = false;     % Flag: true = generate new path, false = load existing

% Path generation or loading
if createPath
   
    % Define paths to dependencies
    % OBS: you must have CasADi installed (https://web.casadi.org/get/)
    addpath(genpath('/home/filipe/Documents/Research/smpc-with-risk-allocation/Example 2/model_in_yop/yop'))
    addpath(genpath('/home/filipe/Documents/casadi-3.7.2'))

    % Generate time-optimal path
    % Calls the external Yop/CasADi function
    path = pathplanner(m1,m2,l,g,ctrl_intervals,ctrl_const);
    
    % Clean up path after generation to avoid conflicts later
    restoredefaultpath;

    % --- Data extraction and resampling ---
    % Construct matrix from planner output
    % Row 1: Independent variable (Time)
    % Row 2-4: States (Velocity, Angle, AngVelocity) 
    % Note: State(1) (Position) is used separately as the interpolation base
    x_path_original = [path.Independent;
        path.State(2,:);
        path.State(3,:);
        path.State(4,:)];

    % Use position (State 1) as the reference for resampling
    time_irr = path.State(1,:);
    
    % Calculate sampling frequency based on total distance/intervals
    fs = ctrl_intervals/time_irr(end);
    Ts = time_irr(end)/ctrl_intervals;
    
    % Resample trajectory to fixed grid (dependenci on Signal Processing Toobox)
    [~,time] = resample(x_path_original(1,:),time_irr,fs);
    
    % Resample States (Time, Velocity, Angle, AngVel) and inputs
    x_path = [resample(x_path_original(1,:),time_irr,fs);
        resample(x_path_original(2,:),time_irr,fs);
        resample(x_path_original(3,:),time_irr,fs);
        resample(x_path_original(4,:),time_irr,fs)];
        
    u_path = resample(path.Input,time_irr,fs);

else
    % Load pre-calculated and resampled path
    load path_stored.mat
end


%%
% addpath(genpath('your_path_to_YALMIP'))
% addpath(genpath('your_path_to_mosek'))
addpath(genpath('/home/filipe/Documents/YALMIP'))
addpath(genpath('/home/filipe/Documents/mosek/11.0/toolbox/r2019b'))
yalmip('clear')

%% Crane model
Ac = [0 1 0 0;
    0 0 m2*g/m1 0;
    0 0 0 1;
    0 0 -(m1+m2)*g/(m1*l) 0];

Bc = [0, 1/m1, 0, -1/(m1*l)]';
Gc = [0, 0, 0, -1]';
Cc = eye(4);
syscc = ss(Ac,[Bc,Gc],Cc,0);

% Discretization
sysdc = c2d(syscc,Ts);
A = sysdc.A;
B = sysdc.B(:,1);
G = sysdc.B(:,2);
C = sysdc.C;

% MPC parameters
nx = size(A,1); % Number of states
nu = size(B,2); % Number of inputs
ny = size(C,1); % Output dimension
nw = size(G,2); % Disturbance dimension

Q = diag([1 1 0 0]); % reference following trolley position and speed only
N = 10; % horizon

% Stochastic parameters
mu = 0.2;
sigma = 0.1;

uppc = 0.15;
lowc = -0.15;

%% Modeling the optimization problem

% Decision variables
u = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % inputs
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1)); % states
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1)); % reference
w = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % disturbance
gammaU = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % Gamma (upper bound)
gammaL = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % Gamma (lower bound)

% Objective and constraint construction
objective = 0;
constraints = [];
gammaSum = 0;

% The controller

for i=1:N
    % System Dynamics
    x{i+1} = A*x{i}+B*u{i}+G*w{i};

    % Sum of the risk allocation
    gammaSum = gammaSum+gammaU{i}+gammaL{i};

    % Constraints (Input, Chance, Uncertainty)
    constraints = [constraints, ...
        -2<=u{i}<=2, ...
        gammaU{i} >= 0, ...
        gammaL{i} >= 0, ...
        probability(x{i+1}(3,1) <= uppc) >= 1-gammaU{i}, ...
        probability(x{i+1}(3,1) >= lowc) >= 1-gammaL{i}, ...
        uncertain(w{i},'normal',mu,sigma)];

    % Objective function
    % Note: Replace w with 0 for nominal cost calculation
    y = replace(C*x{i}-r{i},[w{:}],0);
    objective = objective + y'*Q*y + gammaU{i}'*gammaU{i} + gammaL{i}'*gammaL{i};
end

% Terminal Cost and Final Constraints
constraints = [constraints, 0 <= gammaSum <= 0.15];
y = replace(C*x{N+1}-r{N+1},[w{:}],0);
objective = objective + y'*Q*y;

% Create optimizer object
parameters_in = {x{1},[r{:}]};
solutions_out = {[u{:}],replace([x{:}],[w{:}],0),[gammaU{:}],[gammaL{:}]};

% define the controller
ops = sdpsettings('chance.expcone','yes','solver','mosek','verbose',0);
controller = optimizer(constraints, objective, ops,parameters_in,solutions_out);

%% Execution

% Generate disturbance realization
rng('shuffle');
disturbance = mu+chol(sigma)*randn(1,size(time,2)-1);

% Simulation Initialization
xx = [0;0;0;0];
xhist = xx;
uhist = [];
gammahistU = [];
gammahistL = [];

for k = 1:(size(x_path,2)-1)

    % Construct future reference
    if k<=(size(x_path,2)-N)
        future_r = x_path(:,k:k+N);
    else
        % Reference if near end of path
        future_r = [x_path(:,k:end) repmat(x_path(:,end),1,(N-size(x_path(1:2,k:end),2)+1))];
    end

    % Call controller
    inputs = {xx,future_r};
    [solutions,diagnostics] = controller{inputs};

    if diagnostics == 1
        error('The problem is infeasible');
    end

    % Extract solutions
    U = solutions{1};
    X = solutions{2};
    GammaU = solutions{3};
    GammaL = solutions{4};



    % Store history
    uhist = [uhist, U(1)];
    gammahistU = [gammahistU, GammaU'];
    gammahistL = [gammahistL, GammaL'];

    % State update (apply input to the system)
    xx = A*xx + B*U(1)+G*disturbance(k);
    xhist = [xhist xx];
end

%% Plots
close all

% Color pallete
cor = '#94FF03';
cot = '#2CABFA';
cou = '#9C46FF';

% --- Plot 1: Trajectory Tracking ---
f = figure('Position', [0, 50, 1600, 900]);

% Subplot 1: Trolley Position
subplot(311); hold on
plot(time,x_path(1,:),'color',cor,'LineWidth',1)
plot(time,xhist(1,:),'color',cot,'LineWidth',2)
axis([time(1) time(end) 0 3.1])
grid on
f1 = gca;
f1.XTickLabel = [];
f1.YAxis.FontSize = 20;
f1.GridLineStyle = '-.';
ylabel('$x_c$','Interpreter','latex','FontSize',30)

% Subplot 2: Trolley Velocity
subplot(312); hold on
plot(time,x_path(2,:),'color',cor,'LineWidth',1)
plot(time,xhist(2,:),'color',cot,'LineWidth',2)
axis([time(1) time(end) 0 1.525])
grid on
f2 = gca;
f2.XTickLabel = [];
f2.YAxis.FontSize = 20;
f2.GridLineStyle = '-.';
ylabel('$\dot{x}_c$','Interpreter','latex','FontSize',30)

% Subplot 3: Sway Velocity
subplot(313); hold on
plot(time,xhist(4,:),'color',cot,'LineWidth',2)
grid on
f3 = gca;
f3.YAxis.FontSize = 20;
f3.GridLineStyle = '-.';
f3.XAxis.FontSize = 20;
axis([time(1) time(end) -0.4 0.4205])
xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',30)

% --- Plot 2: Analysis using colormap (Inputs, Angle, Gammas) ---
gidx = 1:N;
cp = clmap(6);
cmap = colormap(cp/255);

b = figure('Position', [0, 50, 1600, 1000]);
tlo = tiledlayout(4,1);

% Tile 1: Control Inputs
nexttile(tlo);hold on
stairs(time(1:end-1),uhist(:),'color',cou,'LineWidth',2)
stairs(time(1:end-1),u_path(1:end-1),'-.k','LineWidth',2)
legend('$u$','$u^r$','fontsize',25,'Interpreter','latex','Location','Eastoutside')
legend boxoff
grid on
axis([time(1) time(end) -2.15 2.15])
j1 = gca;
j1.XTickLabel = [];
j1.YAxis.FontSize = 20;
j1.Box = 'off';
j1.GridLineStyle = '-.';

% Tile 2: Sway Angle
nexttile(tlo); hold on
plot(time,x_path(3,:),'k','LineWidth',2)
plot(time,xhist(3,:),'color',cot,'LineWidth',2)
grid on
axis([time(1) time(end) -0.199999 0.199999])
g1 = gca;
g1.XTickLabel = [];
g1.YAxis.FontSize = 20;
g1.Box = 'off';
g1.GridLineStyle = '-.';
yline(uppc,'k',{'$\theta_{max}$'},'Interpreter','latex','LineWidth',1.5,'FontSize',15, ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','top')
yline(lowc,'k',{'$\theta_{min}$'},'Interpreter','latex','LineWidth',1.5,'FontSize',15, ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','top')
ylabel('$\theta$','Interpreter','latex','FontSize',30)

% Tile 3: Gamma Upper
nexttile(tlo);
imagesc(h(1),time(2:end),gidx,gammahistU(1:N,:));
axis([time(1) time(end) 1 10.5])
set(gca, 'YDir','normal')
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
h1.XTickLabel = [];
ylabel('$i$','Interpreter','latex','FontSize',30)
title('$\gamma^u$','Interpreter','Latex','FontSize',25)

% Tile 4: Gamma Lower
nexttile(tlo);
imagesc(h(2),time(2:end),gidx,gammahistL(1:N,:))
axis([time(1) time(end) 1 10.5])
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
set(gca, 'YDir','normal')
xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$i$','Interpreter','latex','FontSize',30)

% Colorbar adjustments
set(h, 'Colormap', cmap)
cbh = colorbar(h(end));
cbh.Position = [0.93 0.109 0.013 0.375];
cbh.FontSize = 16;
cbh.Label.Interpreter = 'Latex';
cbh.Label.FontSize = 20;
cbh.Label.String = '$\gamma$';
cbh.Label.Rotation = 0;
title('$\gamma^l$','Interpreter','Latex','FontSize',25)
