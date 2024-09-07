% 绘制建筑位置分布图
function PlotBuildingLocations(x, y, type)
    scatter(x, y, 'filled');
    xlabel('横坐标 (m)');
    ylabel('纵坐标 (m)');
    title('建筑位置分布图');
     hold on;
    % 在每个建筑位置上标注建筑编号
    text_offset = 45; % 调整文本偏移量，可以根据需要调整
    for i = 1:length(x)
        text(x(i), y(i) + text_offset, num2str(i), 'HorizontalAlignment', 'center');
    end
end

