% This code is for the reproduction of the paper:
%
% Marques Barbosa, Filipe, and Johan Löfberg. 2025. "Exponential Cone 
% Approach to Joint Chance Constraints in Stochastic Model Predictive 
% Control." International Journal of Control.
% DOI: https://doi.org/10.1080/00207179.2025.2492305

addpath(genpath('your_path_to_YALMIP'))
addpath(genpath('your_path_to_mosek'))
yalmip('clear')

% Model
A = [-0.1 -0.3; 1 -0.5];
B = [2; 0];
G = [1; 0];
C = [1 0];

nx = size(A,1); % number of states
nu = size(B,2); % number of inputs
nw = size(G,2); % dimesion of disturbance
ny = size(C,1); % output dimension

N = 6; % horizon

% The variables
u = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % inputs
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1)); % states
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1)); % reference
w = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % disturbances
gammaU = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % risk of upper bound violation
gammaL = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % risk of lower bound violation

% declare objective and constraints
objective = 0;
constraints = [];

% Stochastic part
gammaSum = 0; % upper bound for Boole's innequality
mu = 0.1; % distribution mean
sigma = 0.2; % distribution covariance
uppc = 3.4; % upper bound
lowc = -3.4; % lower bound

%% The controller
for i=1:N
    
    x{i+1} = A*x{i}+B*u{i}+G*w{i};
    
    gammaSum = gammaSum+gammaU{i}+gammaL{i};

    constraints = [constraints, 
        -2<=u{i}<=2, % input constraints
        gammaU{i} >= 0, % risk cannot be negative
        gammaL{i} >= 0,
        probability(x{i+1}(1,1) <= uppc) >= 1-gammaU{i}, % probability to violate the upper bound constraints
        probability(x{i+1}(1,1) >= lowc) >= 1-gammaL{i}, % probability to violate the lower bound constraints
        uncertain(w{i},'normal',mu,sigma)];

    y = replace(C*x{i}-r{i},[w{:}],0); % cost cannot have stochastic variables

    % cost is quadratic trade-off between risk and performance
    objective = objective + y'*y + gammaU{i}'*gammaU{i} + gammaL{i}'*gammaL{i};
end

constraints = [constraints,
    0<=gammaSum<=0.5]; % bound on joint constraints

y = replace(C*x{N+1}-r{N+1},[w{:}],0);
objective = objective + y'*y;

% Construction the controller
parameters_in = {x{1},[r{:}]};
solutions_out = {[u{:}],replace([x{:}],[w{:}],0),[gammaU{:}],[gammaL{:}]};

ops = sdpsettings('chance.expcone','yes','solver','mosek','verbose',0);
controller = optimizer(constraints, objective,ops,parameters_in,solutions_out);

%% Execution

% some initializations
rng('default');
time = 0:0.0250:6.25; % time vector
disturbance = mu+chol(sigma)*randn(1,size(time,2)-1); % computing the disturbance realization
rt = disturbance;

xx = [0;0]; % initial state
xhist = xx;
uhist = [];
gammahistU = [];
gammahistL = [];
solthist = [];

% For the Monte Carlo simulation
monteCarloSimulations = true;
nsim = 1e7; % Number of realizations

% evaluations points
te1 = 45; 
te2 = 81;
te3 = 171;
te4 = 213;

for k = 1:(size(time,2)-1)
    future_r = 3*sin((k:k+N)/40); % reference to be followed

    inputs = {xx,future_r};
    [solutions,diagnostics,~,~,~,solt] = controller{inputs};

    % extract the solution
    U = solutions{1};
    X = solutions{2};
    GammaU = solutions{3};
    GammaL = solutions{4};
    solthist = [solthist,solt.solvertime];

    if diagnostics == 1
        error('The problem is infeasible');
    end

    % Monte Carlo simulations
    if monteCarloSimulations
        if k==te1 || k==te2 || k==te3 || k==te4
            [individual_probabilities_upper,individual_probabilities_low, joint_probs] = ...
                probabilityMonteCarlo(X(:,2:end),nsim,mu,sigma,GammaU,GammaL,uppc,lowc,N,nx,nw,A,G)
        end
    end

    % Storing the results
    uhist = [uhist, U(1)];
    gammahistU = [gammahistU, GammaU'];
    gammahistL = [gammahistL, GammaL'];

    xx = A*xx + B*U(1)+G*disturbance(k);
    xhist = [xhist xx];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots

reference = 3*sin(time(2:end));
% colors for reference and performed trajectory
cor = '#94FF03';
cot = '#2CABFA';

% colors for upper bound gammas
cog1 = '#F52E12';
cog2 = '#F57119';
cog3 = '#F4E62D';
cog4 = '#F5A715';
cog5 = '#8CF515';
cog6 = '#24F577';


% colors for lower bound gammas
cog7 = '#F53727';
cog8 = '#F57228';
cog9 = '#753D38';
cog10 = '#A06D4F';
cog11 = '#278AF5';
cog12 = '#27F5F4';

% colors for the feasible sets
co0 = [0.0549 0.4941 0.9216];
co2 = [0.9216 0.4549 0.0549];
gray = [0.3 0.3 0.3];

figure('Position', [0, 100, 1600, 900]);
title('Performed trajectory and reference','Interpreter','Latex','FontSize',25)
hold on
plot(time(2:end),reference,'color',cor,'LineWidth',2)
plot(time,xhist(1,:),'color',cot,'LineWidth',3)
xlabel('$t$','Interpreter', 'latex','FontSize',35)
xline([time(te1) time(te2) time(te3) time(te4)],'--k',{'t = 1.10','t = 2.00','t = 4.25','t = 5.30'},...
    'LineWidth',1,'FontSize',14,'LabelHorizontalAlignment','center','LabelOrientation','horizontal')
yline(uppc,'k',{'Upper bound'},'LineWidth',4,'FontSize',20, ...
    'LabelHorizontalAlignment','center','LabelVerticalAlignment','top')
yline(lowc,'k',{'Lower bound'},'LineWidth',4,'FontSize',20, ...
    'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom')
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
legend('$r$','$y$','fontsize',25,'Interpreter','latex','Location','Northeast')
legend boxoff
axis([0 6.25 -5 5])


figure('Position', [0, 100, 1600, 900]);
nexttile;
title('Confidence parameters','Interpreter','Latex','FontSize',25)
hold on
plot(time(2:end),gammahistU(1,:),'color',cog1,'LineWidth',2)
plot(time(2:end),gammahistU(2,:),'color',cog2,'LineWidth',2)
plot(time(2:end),gammahistU(3,:),'color',cog3,'LineWidth',2)
plot(time(2:end),gammahistU(4,:),'color',cog4,'LineWidth',2)
plot(time(2:end),gammahistU(5,:),'color',cog5,'LineWidth',2)
plot(time(2:end),gammahistU(6,:),'color',cog6,'LineWidth',2)
plot(time(2:end),gammahistL(1,:),'-.','color',cog1,'LineWidth',2)
plot(time(2:end),gammahistL(2,:),'-.','color',cog2,'LineWidth',2)
plot(time(2:end),gammahistL(3,:),'-.','color',cog3,'LineWidth',2)
plot(time(2:end),gammahistL(4,:),'-.','color',cog4,'LineWidth',2)
plot(time(2:end),gammahistL(5,:),'-.','color',cog5,'LineWidth',2)
plot(time(2:end),gammahistL(6,:),'-.','color',cog6,'LineWidth',2)
xline([time(te1) time(te2) time(te3) time(te4)],'--k',{'t = 1.10','t = 2.00','t = 4.25','t = 5.30'},...
    'LineWidth',1,'FontSize',14,'LabelHorizontalAlignment','center','LabelOrientation','horizontal')
xlabel('$t$','Interpreter', 'latex','FontSize',35)
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
leg = legend;
legend('$\gamma^{u}_{k+1}$','$\gamma^{u}_{k+2}$','$\gamma^{u}_{k+3}$','$\gamma^{u}_{k+4}$','$\gamma^{u}_{k+5}$','$\gamma^{u}_{k+6}$', ...
    '$\gamma^{l}_{k+1}$','$\gamma^{l}_{k+2}$','$\gamma^{l}_{k+3}$','$\gamma^{l}_{k+4}$','$\gamma^{l}_{k+5}$','$\gamma^{l}_{k+6}$', ...
    'fontsize',25,'Interpreter','latex','Location','north','Orientation','Vertical')
leg.NumColumns = 2;
legend boxoff
axis([0 6.25 0 0.2])


%% Plot feasible set
figure('Position', [200, 100, 1600, 900]);
hold on
title('Feasible sets','Interpreter','Latex','FontSize',25)
plot(constraints,x{1},co0,100,sdpsettings('chance.expcone','yes','plot.edgecolor','none','plot.shade',1))
plot([constraints,[gammaL{:}]==0.5/(2*N),[gammaU{:}]==0.5/(2*N)],x{1},co2,100,...
    sdpsettings('chance.expcone','yes','plot.edgecolor','none','plot.shade',1))
axis([-32 32 -30 30])
h = gca;
h.XColor = gray;
h.YColor = gray;
h.XAxis.FontSize = 20;
h.YAxis.FontSize = 20;
xlabel('$x_1$','Interpreter','Latex','FontSize',35)
ylabel('$x_2$','Interpreter','Latex','FontSize',35)
legend('Proposed approach','Bonferonni correction','Interpreter','latex','FontSize',20,'Location','Northeast');
legend boxoff;
