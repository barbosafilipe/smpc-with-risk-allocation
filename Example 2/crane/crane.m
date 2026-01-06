%% Obtaining the reference path
ctrl_intervals = 100;
m1 = 1.2;
m2 = 0.6;
l = 1.5;
g = 9.8067;
ctrl_const = [-2,2];
% Uncomment bellow to generate a new path 
% path = pathplanner(m1,m2,l,g,ctrl_intervals,ctrl_const);
% restoredefaultpath;

% use a predefined path (they gotta be always the same)
% load path_stored.mat

x_path_original = [path.Independent;
    path.State(2,:);
    path.State(3,:);
    path.State(4,:)];
time_irr = path.State(1,:);
fs = ctrl_intervals/time_irr(end);
Ts = time_irr(end)/ctrl_intervals;

[~,time] = resample(x_path_original(1,:),time_irr,fs);

x_path = [resample(x_path_original(1,:),time_irr,fs);
    resample(x_path_original(2,:),time_irr,fs);
    resample(x_path_original(3,:),time_irr,fs);
    resample(x_path_original(4,:),time_irr,fs)];

u_path = resample(path.Input,time_irr,fs);

%%
addpath(genpath('C:\Users\filma90\Documents\MATLAB\YALMIP-2023\YALMIP'))
addpath(genpath('C:\Program Files\Mosek\10.1\toolbox\r2017a'))
yalmip('clear')

%% Crane model
Ac = [0 1 0 0; 0 0 m2*g/m1 0; 0 0 0 1; 0 0 -(m1+m2)*g/(m1*l) 0];
Bc = [0, 1/m1, 0, -1/(m1*l)]';
Gc = [0, 0, 0, -1]';
Cc = eye(4);
syscc = ss(Ac,[Bc,Gc],Cc,0);
sysdc = c2d(syscc,Ts);
A = sysdc.A;
B = sysdc.B(:,1);
G = sysdc.B(:,2);
C = sysdc.C;

nx = 4; % Number of states
nu = 1; % Number of inputs
ny = 4; % output dimension

% MPC data
Q = diag([1 1 0 0]); % reference following only the trolley position and speed
N = 10; % horizon


% Obtaining a linear map
x_base = sdpvar(nx,1);
u_base = sdpvar(nu,N);
w_base = sdpvar(1,N);
X_base = [];

for k = 1:N
    x_base = A*x_base + B*u_base(:,k) + G*w_base(:,k);
    X_base = [X_base; x_base];
end
lmap = full(getbase(X_base));
lmap = lmap(:,2:end);

Amap = lmap(:,1:nx);
Bmap = lmap(:,(nx+1):(nx+N));
Gmap = lmap(:,(nx+1+N):(nx+2*N));

u = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % inputs
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1)); % states
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1)); % reference
w = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % disturbance
gammaU = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % gamma for upper bound
gammaL = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % gamma for lower bound

% declare objective and constraints
objective = 0;
constraints = [];

%% Stochastic part
mu = 0.2;
sigma = 0.1;
gammaSum = 0;

uppc = 0.15;
lowc = -0.15;

%% The controller

for i=1:N
    x{i+1} = A*x{i}+B*u{i}+G*w{i};
    gammaSum = gammaSum+gammaU{i}+gammaL{i};
    constraints = [constraints, -2<=u{i}<=2,
                   gammaU{i} >= 0,
                   gammaL{i} >= 0,
                   probability(x{i+1}(3,1) <= uppc) >= 1-gammaU{i},
                   probability(x{i+1}(3,1) >= lowc) >= 1-gammaL{i},
                   uncertain(w{i},'normal',mu,sigma)];
    
    y = replace(C*x{i}-r{i},[w{:}],0);
    objective = objective + y'*Q*y + gammaU{i}'*gammaU{i} + gammaL{i}'*gammaL{i};
end
constraints = [constraints, 0<=gammaSum<=0.15];
y = replace(C*x{N+1}-r{N+1},[w{:}],0);
objective = objective + y'*Q*y;

parameters_in = {x{1},[r{:}]};
solutions_out = {[u{:}],replace([x{:}],[w{:}],0),[gammaU{:}],[gammaL{:}]};

% define the controller
ops = sdpsettings('chance.expcone','yes','solver','mosek','verbose',0);
controller = optimizer(constraints, objective, ops,parameters_in,solutions_out);

%% Execution
rng('shuffle');
disturbance = mu+chol(sigma)*randn(1,size(time,2)-1);

figure(1)
hold on
plot(time(2:end),disturbance,'.')
axis([0 3.5 -3 3])
%%
xx = [0;0;0;0];
xhist = xx;
uhist = [];
gammahistU = [];
gammahistL = [];
nsim = 1e7;
te1 = 78;

for k = 1:(size(x_path,2)-1)
    
    if k<=(size(x_path,2)-N)
        future_r = x_path(:,k:k+N);
    else
        future_r = [x_path(:,k:end) repmat(x_path(:,end),1,(N-size(x_path(1:2,k:end),2)+1))];
    end

    inputs = {xx,future_r};
    [solutions,diagnostics] = controller{inputs};
    U = solutions{1};
    X = solutions{2};
    GammaU = solutions{3};
    GammaL = solutions{4};

    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    % chance constraints
    uhist = [uhist, U(1)];
    gammahistU = [gammahistU, GammaU'];
    gammahistL = [gammahistL, GammaL'];

    xx = A*xx + B*U(1)+G*disturbance(k);
    xhist = [xhist xx];
end

%% Plots
close all

% colors for reference and performed trajectory
cor = '#94FF03';
cot = '#2CABFA';
cou = '#9C46FF';

f = figure('Position', [0, 50, 1600, 900]);
subplot(311); hold on
plot(time,x_path(1,:),'color',cor,'LineWidth',1)
plot(time,xhist(1,:),'color',cot,'LineWidth',2)
axis([time(1) time(end) 0 3.1])
f1 = gca;
f1.XTickLabel = [];
f1.YAxis.FontSize = 20;
f1.GridLineStyle = '-.';
grid on
ylabel('$x_c$','Interpreter','latex','FontSize',30)
subplot(312); hold on
plot(time,x_path(2,:),'color',cor,'LineWidth',1)
plot(time,xhist(2,:),'color',cot,'LineWidth',2)
axis([time(1) time(end) 0 1.525])
f2 = gca;
f2.XTickLabel = [];
f2.YAxis.FontSize = 20;
f2.GridLineStyle = '-.';
grid on
ylabel('$\dot{x}_c$','Interpreter','latex','FontSize',30)
subplot(313); hold on
plot(time,xhist(4,:),'color',cot,'LineWidth',2)
f3 = gca;
f3.YAxis.FontSize = 20;
f3.GridLineStyle = '-.';
f3.XAxis.FontSize = 20;
axis([time(1) time(end) -0.4 0.4205])
grid on
xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',30)
drawnow
% exportgraphics(f,'C:\Users\filma90\Documents\Texts\writing\Journal\Chance_constraints\fig\traj-crane.eps')

%% Color map plot


gidx = 1:N;
cp = clmap(6);
cmap = colormap(cp/255);

b = figure('Position', [0, 50, 1600, 1000]);
tlo = tiledlayout(4,1);
j(1) = nexttile(tlo);
hold on
stairs(time(1:end-1),uhist(:),'color',cou,'LineWidth',2)
% stairs(time(1:end-1),u_path(1:end-1),'-.k','LineWidth',2)
legend('$u$','$u^r$','fontsize',25,'Interpreter','latex','Location','Eastoutside')
legend boxoff
grid on
axis([time(1) time(end) -2.15 2.15])
j1 = gca;
j1.XTickLabel = [];
j1.YAxis.FontSize = 20;
j1.Box = 'off';
j1.GridLineStyle = '-.';
g(1) = nexttile(tlo);
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
h(1) = nexttile(tlo);
imagesc(h(1),time(2:end),gidx,gammahistU(1:N,:));
axis([time(1) time(end) 1 10.5])
% xline([time(te1)],'--k','LineWidth',1)
set(gca, 'YDir','normal')
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
h1.XTickLabel = [];
ylabel('$i$','Interpreter','latex','FontSize',30)
title('$\gamma^u$','Interpreter','Latex','FontSize',25)
h(2) = nexttile(tlo); 
imagesc(h(2),time(2:end),gidx,gammahistL(1:N,:))
axis([time(1) time(end) 1 10.5])
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
set(gca, 'YDir','normal')
xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$i$','Interpreter','latex','FontSize',30)
set(h, 'Colormap', cmap)
cbh = colorbar(h(end));
% cbh.Layout.Tile = 'east';
cbh.Position = [0.93 0.109 0.013 0.375];
cbh.FontSize = 16;
cbh.Label.Interpreter = 'Latex';
cbh.Label.FontSize = 20;
cbh.Label.String = '$\gamma$';
cbh.Label.Rotation = 0;
title('$\gamma^l$','Interpreter','Latex','FontSize',25)
% exportgraphics(b,'C:\Users\filma90\Documents\Texts\writing\Journal\Chance_constraints\fig\color-angle.eps')

%%
gidx = 1:N;
cp = clmap(6);
cmap = colormap(cp/255);

b = figure('Position', [0, 50, 1600, 1000]);
tlo = tiledlayout(3,1);
g(1) = nexttile(tlo);
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
h(1) = nexttile(tlo);
imagesc(h(1),time(2:end),gidx,gammahistU(1:N,:));
axis([time(1) time(end) 1 10.5])
% xline([time(te1)],'--k','LineWidth',1)
set(gca, 'YDir','normal')
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
h1.XTickLabel = [];
ylabel('$i$','Interpreter','latex','FontSize',30)
title('$\gamma^u$','Interpreter','Latex','FontSize',25)
h(2) = nexttile(tlo); 
imagesc(h(2),time(2:end),gidx,gammahistL(1:N,:))
axis([time(1) time(end) 1 10.5])
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
set(gca, 'YDir','normal')
xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$i$','Interpreter','latex','FontSize',30)
set(h, 'Colormap', cmap)
cbh = colorbar(h(end));
% cbh.Layout.Tile = 'east';
cbh.Position = [0.93 0.109 0.013 0.522];
cbh.FontSize = 16;
cbh.Label.Interpreter = 'Latex';
cbh.Label.FontSize = 20;
cbh.Label.String = '$\gamma$';
cbh.Label.Rotation = 0;
title('$\gamma^l$','Interpreter','Latex','FontSize',25)
