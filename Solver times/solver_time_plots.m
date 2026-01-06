%% Solver times
if 0
    load('time_ipopt.mat')
    load('time_mosek.mat')
    load('time_fmincon.mat')

    figure('Position', [0, 100, 1600, 900]);
    hold on
    plot(ipopt_time(1,:),'color',[0,0,1],'LineWidth',3)
    plot(mosek_time(1,:),'color',[0,1,0],'LineWidth',3)
    % plot(fmincon_time(1,:),'color',[1,0,0],'LineWidth',3)
    xlabel('$k$','Interpreter', 'latex','FontSize',35)
    h1 = gca;
    h1.XAxis.FontSize = 20;
    h1.YAxis.FontSize = 20;
    legend('ipopt','mosek','fmincon','fontsize',25,'Interpreter','latex','Location','Northeast')
    legend boxoff
end

%% 
load('Stimes_ipopt.mat')
% load('Stimes_mosek.mat')

solverTimes_avg = [];
solverTimes_std = [];

for i = 1:size(solver_times,2)
    solverTimes_avg(i) = mean(solver_times{i});
    solverTimes_std(i) = std(solver_times{i});
end