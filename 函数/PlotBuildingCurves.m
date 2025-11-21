function PlotBuildingCurves(load_curve, solar_curve)
%     % 绘制负荷曲线
t = linspace(0, 24, 48); % 生成横坐标，48个点均匀分布在0到24小时
    plot(t,load_curve,'LineWidth', 1.5);
    hold on;
%     % 绘制光伏曲线
    plot(t,solar_curve,'LineWidth', 1.5);
    hold on;
    net_load_curve=load_curve-solar_curve;
    plot(t,net_load_curve,'LineWidth', 1.5);
    % 设置 x 轴范围和刻度
xlim([0, 24]); % 确保横坐标范围
xticks(0:4:24); % 强制显示 0, 4, 8, 12, 16, 20, 24
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'}); % 设置刻度标签
    xlabel('时刻');
    ylabel('夏季典型日功率（kW）');
    title('建筑');
    legend('负荷功率曲线','光伏发电功率曲线','净负荷曲线')
    hold off;
end


