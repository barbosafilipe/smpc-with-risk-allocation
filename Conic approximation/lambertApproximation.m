% This code is for obtaining the conic convex approximation used in the paper:
%
% Marques Barbosa, Filipe, and Johan Löfberg. 2025. "Exponential Cone 
% Approach to Joint Chance Constraints in Stochastic Model Predictive 
% Control." International Journal of Control.
% DOI: https://doi.org/10.1080/00207179.2025.2492305

addpath(genpath('your_path_to_YALMIP'))
yalmip('clear')

%% Using the Lambert-W function as an exponential cone approximation to the probit function

x0 = 0.0001; x = x0:0.0001:0.5; % Define the region where we want to approxmate

% Initial guesses of the parameters
a0 = 0.70;
b0 = 600;
c0 = -0.75;
k0 = a0*lambertw(b0*x(1))+norminv(1-x(1));

sdpvar a b c k
z = sdpvar(length(x),1);
warmstart([a b c k],[a0 b0 c0 k0]);

obj = 0;
con = -a*lambertw(b*x(1))+k >= norminv(1-x(1));

% The curve fitting problem is formulated as a non-negative least-squares regression
for i = 1:100:length(x)
    fi = -a*lambertw(b*x(i)) + k + c*x(i);
    yi = norminv(1-x(i));
    con = [con, z(i) == (fi-yi), z(i) >= 0];
    obj = obj + (1 + 100*(i==1))*(z(i))^2;
end
optimize(con,obj,sdpsettings('usex0',1,'solver','fmincon'));

%%
a = value(a);
b=value(b);
k=value(k);
c = value(c);

% The values used on the paper:
% a = 0.499492956059166;
% b = 8.082867432374761e+03;
% k = 3.965651977413067;
% c = -1.475743096725997;

x0 = 0.00001;
x = x0:0.0001:1;
lamb=@(x)(-a*lambertw(b*x)+k+c*x);

%% Plots

co1 = [0.8500 0.3250 0.0980];
co2 = '#F09A11';
co3 = '#188AF0';

%% Plot the probit function and the obtained approximations
figure('Position', [200, 100, 1600, 900]);hold on;
hold on
plot(x,norminv(1-x),'color',co2,'LineWidth',4)
plot(x,lamb(x),':','color',co3,'LineWidth',5)
axis([-0.01 1.01 -4.5 4])
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
xlabel('$\gamma$','Interpreter', 'latex','FontSize',55)
legend('$\Phi^{-1}(1-\gamma)$','$\Psi(\gamma)$','Interpreter','latex','FontSize',35)
legend boxoff;

%% Plot the logarithm error
figure('Position', [200, 100, 1600, 900]);hold on;
hold on
plot(x,log(lamb(x)-norminv(1-x)),'color',co1,'LineWidth',4)
axis([-0.01 0.5 -22 0])
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
xlabel('$\gamma$','Interpreter', 'latex','FontSize',55)
legend('$\log e(\gamma)$','Interpreter','latex','FontSize',35,'Location','NorthEast')
legend boxoff;