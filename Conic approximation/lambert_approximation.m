addpath(genpath('C:\Users\filma90\Documents\MATLAB\YALMIP-2023\YALMIP'))
addpath(genpath('C:\Program Files\Mosek\10.1\toolbox\r2017a'))
yalmip('clear')

%% Using Lambert

x0 = 0.0001;
x = x0:0.0001:0.5;

a0 = 0.70;
b0 = 600;
c0 = -0.75;
k0 = a0*lambertw(b0*x(1))+norminv(1-x(1));

sdpvar a b c k
z = sdpvar(length(x),1);
warmstart([a b c k],[a0 b0 c0 k0]);

obj = 0;
Con = -a*lambertw(b*x(1))+k >= norminv(1-x(1));

for i = 1:100:length(x)
    fi = -a*lambertw(b*x(i)) + k + c*x(i);
    yi = norminv(1-x(i));    
    Con = [Con, z(i) == (fi-yi), z(i) >= 0];
    obj = obj + (1 + 10000*(i==1))*(z(i))^2;
    % assign(z(i), value(fi)-yi);
end
optimize(Con,obj,sdpsettings('usex0',1,'solver','ipopt'));

%%
% a = value(a);
% b=value(b);
% k=value(k);
% c = value(c);

a = 0.499492956059166;
b = 8.082867432374761e+03;
k = 3.965651977413067;
c = -1.475743096725997;

x0 = 0.00001;
x = x0:0.0001:1;
lamb=@(x)(-a*lambertw(b*x)+k+c*x);
xn = x0:0.0001:1;
%% Plots

red = [0.8941 0.0235 0.0235];
gray = [0.5 0.5 0.5];
yellow = [0.2667,0.4471,0.7098];
purple = [0.8275,0.2627,0.3059];
darkGreen = [0.1647,0.6706,0.3804];
blue = [0 0 1];
cyan = [0 1 1];
co1 = [0.9290 0.6940 0.1250];
co2 = [0.8500 0.3250 0.0980];
co3 = '#F09A11';
co4 = '#188AF0';

%%
figure('Position', [200, 100, 1600, 900]);hold on;
hold on
plot(x,norminv(1-x),'color',co3,'LineWidth',4)
plot(x,lamb(x),':','color',co4,'LineWidth',5)
axis([-0.01 1.01 -4.5 4])
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
xlabel('$\gamma$','Interpreter', 'latex','FontSize',55)
legend('$\Phi^{-1}(1-\gamma)$','$\Psi(\gamma)$','Interpreter','latex','FontSize',35)
legend boxoff;
% drawnow;
% exportgraphics(gca,'C:\Users\filma90\Documents\Texts\writing\Journal\IJC\fig\approximation.eps')
%%
figure('Position', [200, 100, 1600, 900]);hold on;
hold on
plot(x,log(lamb(x)-norminv(1-x)),'color',co2,'LineWidth',4)
axis([-0.01 0.5 -22 0])
h2 = gca;
h2.XAxis.FontSize = 20;
h2.YAxis.FontSize = 20;
xlabel('$\gamma$','Interpreter', 'latex','FontSize',55)
legend('$\log e(\gamma)$','Interpreter','latex','FontSize',35,'Location','NorthEast')
legend boxoff;
% drawnow;
% exportgraphics(gca,'C:\Users\filma90\Documents\Texts\writing\Journal\IJC\fig\error.eps')


%%

figure('Position', [200, 100, 1600, 900]);hold on;
hold on
plot(xn,norminv(1-xn),'color',co3,'LineWidth',4)
plot(x,lamb(x),':','color',co4,'LineWidth',5)
axis([-0.01 3.01 -4.5 4])
grid on
h1 = gca;
h1.XAxis.FontSize = 20;
h1.YAxis.FontSize = 20;
h1.GridLineStyle = "-.";
xlabel('$\zeta$','Interpreter', 'latex','FontSize',35)
legend('$\Phi^{-1}(1-\zeta)$','$\Psi(\zeta)$','Interpreter','latex','FontSize',35)
legend boxoff;
% drawnow;
% exportgraphics(gca,'C:\Users\filma90\Documents\Texts\writing\Journal\TCST\fig\approximation-extended.eps')
