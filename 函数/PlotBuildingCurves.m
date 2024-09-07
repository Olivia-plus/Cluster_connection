function PlotBuildingCurves(load_curve, solar_curve)
%     % 绘制负荷曲线
%     subplot(3, 1, 1);
%     plot(load_curve, 'color','#0072BD', 'LineWidth', 1.5);
%     xlabel('时间');
%     ylabel('功率（kW）');
%     title('负荷曲线');
% 
%     % 绘制光伏曲线
%     subplot(3, 1, 2);
%     plot(solar_curve, 'color', '#D95319', 'LineWidth', 1.5);
%     xlabel('时间');
%     ylabel('功率（kW）');
%     title('光伏曲线');

%     subplot(3, 1, 3);
    net_load_curve=load_curve-solar_curve;
    plot(net_load_curve,  'color','#7E2F8E', 'LineWidth', 1.5);
    xlabel('时间');
    ylabel('功率（kW）');
    title('净负荷曲线');

%     % 调整图的布局
%     sgtitle('建筑曲线图');
%     set(gcf, 'Position', [100, 100, 800, 800]);
end


