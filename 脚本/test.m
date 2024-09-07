% n = 3; % 矩阵的大小
% total_matrices = 2^(n*(n-1)/2); % 总可能的矩阵数量
% 
% % 循环遍历所有可能的情况
% for k = 0:total_matrices-1
%     % 将索引 k 转换成二进制，表示当前矩阵的状态
%     binary = dec2bin(k, n*(n-1)/2);
%     idx = 1;
%     
%     % 生成对称矩阵
%     current_matrix = zeros(n);
%     for i = 1:n
%         for j = i+1:n
%             % 检查是否在对角线上，是则置0
%             if i == j
%                 current_matrix(i, j) = 0;
%             else
%                 % 根据二进制值设置当前位置的值
%                 current_matrix(i, j) = str2double(binary(idx));
%                 current_matrix(j, i) = str2double(binary(idx));
%                 idx = idx + 1;
%             end
%         end
%     end
%     
%     % 对当前矩阵进行处理，可以在这里进行你需要的操作
%     disp(current_matrix);
%     % 在这里添加你的处理代码，对current_matrix进行操作
% end

%% 定义一个行向量
% vector = [1, 2, 3, 1, 1, 2, 3, 3, 2];
% 
% % 使用 unique() 函数找出唯一元素和每个元素在原始向量中的位置
% [unique_elements, ~, idx] = unique(vector);
% 
% % 遍历每个唯一元素，并找出原始向量中具有相同元素的项
% for i = 1:numel(unique_elements)
%     element = unique_elements(i);
%     indices = find(idx == i);
%     % 使用 fprintf() 函数输出字符串和数组
%     fprintf('数字 %d 的索引为: %s\n', element, num2str(indices));
% end

% % 定义一个行向量
% vector = [1, 2, 3, 1, 1, 2, 3, 3, 2];
% 
% % 使用 unique() 函数找出唯一元素和每个元素在原始向量中的位置
% [unique_elements, ~, idx] = unique(vector);
% 
% % 创建一个 cell 数组用于存储每个数字的索引
% index_cell = cell(1, numel(unique_elements));
% 
% % 遍历每个唯一元素，并找出原始向量中具有相同元素的项
% for i = 1:numel(unique_elements)
%     element = unique_elements(i);
%     % 找出当前元素的索引
%     indices = find(vector == element);
%     % 将索引存储到 cell 数组中
%     index_cell{i} = indices;
% end
% 
% % 创建一个新的向量，存储数值一样的编号
% new_vector = zeros(size(vector));
% for i = 1:numel(index_cell)
%     indices = index_cell{i};
%     new_vector(indices) = i;
% end
% 
% % 输出新的向量和原始向量中数值一样的编号的映射
% disp('新的向量:');
% disp(new_vector);
% disp('原始向量中数值一样的编号的映射:');
% disp(cell2mat(index_cell));

% % 创建一个元胞数组，每个元素都是一个行向量
% cell_array = {1:3, 4:7, 8:10};
% 
% % 使用 cellfun 函数对每个元素应用 length 函数，得到每个元素的大小
% sizes = cellfun(@length, cell_array);
% 
% % 显示每个元素的大小
% disp('每个元素的大小：');
% disp(sizes);


% % 示例数据
% num_buildings = 20;
% hours_per_day = 24;
% 
% % 假设 daily_load_data 是一个含有 20 个元素的元胞数组，每个元素是一个 1x24 的行向量
% daily_load_data = cell(1, num_buildings);
% for i = 1:num_buildings
%     daily_load_data{i} = rand(1, hours_per_day); % 随机生成示例数据
% end
% 
% % 创建小时和建筑索引网格
% [hour, building] = meshgrid(1:hours_per_day, 1:num_buildings);
% 
% % 将每个建筑的负荷数据放在网格中的对应位置
% load_grid = zeros(num_buildings, hours_per_day);
% for i = 1:num_buildings
%     load_grid(i, :) = daily_load_data{i};
% end
% 
% % 绘制三维图
% figure;
% surf(building, hour, load_grid);
% xlabel('建筑');
% ylabel('小时');
% zlabel('负荷');
% title('每个建筑的日负荷数据');
% colorbar; % 添加颜色条

% % 示例数据
% num_buildings = 20;
% points_per_hour = 48;
% 
% % 假设 daily_load_data 是一个含有 20 个元素的元胞数组，每个元素是一个 1x48 的行向量
% daily_load_data = cell(1, num_buildings);
% for i = 1:num_buildings
%     daily_load_data{i} = rand(1, points_per_hour); % 随机生成示例数据
% end
% 
% % 创建时间轴标签
% time_labels = cell(1, points_per_hour);
% for i = 1:points_per_hour
%     hour = floor((i-1)/2);
%     minute = rem(i-1, 2) * 30;
%     time_labels{i} = sprintf('%02d:%02d', hour, minute);
% end
% 
% % 创建小时和建筑索引网格
% [point, building] = meshgrid(1:points_per_hour, 1:num_buildings);
% 
% % 将每个建筑的负荷数据放在网格中的对应位置
% load_grid = zeros(num_buildings, points_per_hour);
% for i = 1:num_buildings
%     load_grid(i, :) = daily_load_data{i};
% end
% 
% % 绘制三维图
% figure;
% surf(building, point, load_grid);
% xlabel('建筑');
% ylabel('时间');
% zlabel('负荷');
% title('每个建筑的日负荷数据');
% 
% % 设置时间轴标签
% xticks(1:num_buildings);
% xticklabels(1:num_buildings);
% yticks(1:6:points_per_hour);
% yticklabels(time_labels(1:6:points_per_hour));

% 第一种情况的数据
data1 = [654.2546, 2428.6906, 3082.9452, 6468.9089, 51.1517, 21.3878, 204.8654, 277.4049];

% 第二种情况的数据
data2 = [791.3076, 251.5046, 1042.8122, 6365.8616, 21.9631, 15.4273, 166.1144, 203.5048];

% 将数据转换为列向量
data1 = data1';
data2 = data2';

% 绘制柱状图
bar([data1, data2]);
legend('Scenario 1', 'Scenario 2');
xlabel('Parameters');
ylabel('Values');
title('Comparison of Scenarios');

%%
% 储能容量配置数据
capacity_A = [654.2546, 2428.6906];
capacity_B = [791.3076, 251.5046];

% 绘制柱状图
figure;
bar([capacity_A; capacity_B], 'stacked');
xlabel('两种情况');
ylabel('储能容量 (kWh)');
title('建筑储能容量配置');
legend('Building A', 'Building B');

% 添加标注
text(1, capacity_A(1), num2str(capacity_A(1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(1, sum(capacity_A), num2str(sum(capacity_A)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

text(2, capacity_A(2), num2str(capacity_A(2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(2, sum(capacity_A)+capacity_B(1), num2str(capacity_B(1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(2, sum(capacity_A)+capacity_B(1), num2str(capacity_B(1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(2, sum(capacity_A)+sum(capacity_B), num2str(sum(capacity_B)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

set(gca, 'XTickLabel', {'不互联', '互联'});


%%
% 数据
invest_cost = [51.1517, 21.9631];
maintenance_cost = [21.3878, 15.4273];
grid_cost = [204.8654, 166.1144];
total_cost = [277.4049, 203.5048];

% 绘制堆叠柱状图
figure;
bar([invest_cost; maintenance_cost; grid_cost; total_cost]', 'stacked');
xlabel('Scenarios');
ylabel('Cost');
title('Comparison of Costs');
legend('Investment', 'Maintenance', 'Grid Purchase', 'Total');

% 添加标签
xticks([1 2]);
xticklabels({'Scenario 1', 'Scenario 2'});

% 设置颜色
colors = [0.85 0.33 0.1; 0.93 0.69 0.13; 0.47 0.67 0.19; 0.3 0.75 0.93];
colormap(colors);

%%
% 数据
invest_cost = [51.1517, 21.9631];
maintenance_cost = [21.3878, 15.4273];
grid_cost = [204.8654, 166.1144];
total_cost = [277.4049, 203.5048];
grid_power = [6468.9089, 6365.8616]; % 电网购电总功率数据

% 绘制分组柱状图
figure;
bar([invest_cost; maintenance_cost; grid_cost; total_cost; grid_power]');
xlabel('Cost Components');
ylabel('Value(万元)');
title('成本及综合费用对比');

% 截断函数
truncAxis('Y',[1000,6000])

% 添加标签
legend('Investment ', 'Investment ', ...
       'Maintenance ', 'Maintenance ', ...
       'Grid Purchase ', 'Grid Purchase ', ...
       'Total ', 'Total ', ...
       'Grid Power', 'Grid Power', 'Location', 'best');

% 调整标签位置
xticks(1:5);
xticklabels({'Investment', 'Maintenance', 'Grid Purchase', 'Total', 'Grid Power'});
% 调整字体
set(gca, 'FontSize', 12);
grid on

%%
% 数据
invest_cost = [51.1517, 21.9631];
maintenance_cost = [21.3878, 15.4273];
grid_cost = [204.8654, 166.1144];
total_cost = [277.4049, 203.5048];
grid_power = [6468.9089, 6365.8616]; % 电网购电总功率数据

% 合并数据
data = [invest_cost; maintenance_cost; grid_cost; total_cost; grid_power];

% 绘制分组柱状图
figure;
bar(data);
xlabel('Cost Components');
ylabel('Value(万元)');
title('成本及综合费用对比');

% 截断函数
truncAxis('Y',[400,6000])

% 添加标签
legend('Investment (Scenario 1)', 'Investment (Scenario 2)', ...
       'Maintenance (Scenario 1)', 'Maintenance (Scenario 2)', ...
       'Grid Purchase (Scenario 1)', 'Grid Purchase (Scenario 2)', ...
       'Total (Scenario 1)', 'Total (Scenario 2)', ...
       'Grid Power (Scenario 1)', 'Grid Power (Scenario 2)', 'Location', 'best');

% 调整标签位置
xticks(1:5);
xticklabels({'Investment', 'Maintenance', 'Grid Purchase', 'Total', 'Grid Power'});
xtickangle(45); % 旋转标签，使其更易读

% 调整字体
set(gca, 'FontSize', 12);
grid on;

%% 2个建筑
% 数据
invest_cost = [51.1517, 21.9631];
maintenance_cost = [21.3878, 15.4273];
grid_cost = [204.8654, 166.1144];
total_cost = [277.4049, 203.5048];
grid_power = [6468.9089, 6365.8616]; % 电网购电总功率数据

% 合并数据
data = [invest_cost; maintenance_cost; grid_cost; total_cost; grid_power];

% 绘制分组柱状图
figure;
bar(data);
xlabel('Cost Components');
ylabel('Value(万元)');
title('成本及综合费用对比');

% 截断函数
truncAxis('Y',[1000,6000])

% 添加图例

% 调整标签位置
xticks(1:5);
xticklabels({'设备投资成本', '设备维护成本', '电网购电成本', '系统综合成本', '电网购电总功率'});
xtickangle(45); % 旋转标签，使其更易读

% 调整字体
set(gca, 'FontSize', 12);
grid on;

%% 3个建筑
capacity_A = [654.2546, 3892.8539];
capacity_B = [1758.6845, 1222.6278];
capacity_C = [4418.9177, 0];

% 合并数据
data = [capacity_A; capacity_B; capacity_C];

% 绘制柱状图
figure;
bar(data', 'stacked');
xlabel('两种情况');
ylabel('储能容量 (kWh)');
title('建筑储能容量配置');
legend('Building A', 'Building B', 'Building C', 'Location', 'best');

% 添加标注
for i = 1:numel(data)
    text(i, data(i), num2str(data(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% 设置 x 轴刻度标签
set(gca, 'XTickLabel', {'不互联', '互联'});

%% 2个建筑
% 数据
invest_cost = [346.9576, 287.3049];
maintenance_cost = [183.8749, 126.8468];
grid_cost = [236.0179, 144.0865];
total_cost = [766.8504, 558.2382];
grid_power = [10054.7797,6825.6845]; % 电网购电总功率数据

% 合并数据
data = [invest_cost; maintenance_cost; grid_cost; total_cost; grid_power];

% 绘制分组柱状图
figure;
bar(data);
xlabel('Cost Components');
ylabel('Value(万元)');
title('成本及综合费用对比');

% 截断函数
truncAxis('Y',[1000,6500])

% 添加图例

% 调整标签位置
xticks(1:5);
xticklabels({'设备投资成本', '设备维护成本', '电网购电成本', '系统综合成本', '电网购电总功率'});
xtickangle(45); % 旋转标签，使其更易读

% 调整字体
set(gca, 'FontSize', 12);
grid on;

%% 
[x, y, type,load_curve,pv_curve]=GenerateBuildingInfo(4)

%%

~any(ismember(1:5, particles(i, :)))

%% 主程序 PSO
clear
close all
clc
SearchAgents_no = 30 ; %　种群规模
dim = 10 ; % 粒子维度，函数的变量个数
Max_iter = 1000 ; %　迭代次数
ub = 5 ;
lb = -5 ;
c1 = 1.5 ; %　学习因子1
c2 = 1.5 ; %　学习因子2
w = 0.8 ; %　惯性权重
vmax = 3 ; %　最大飞行速度
pos = lb + rand(SearchAgents_no,dim).*(ub-lb) ; % 初始化粒子群的位置
v = - vmax +2*vmax* rand(SearchAgents_no,dim) ; % 初始化粒子群的速度【rand()生成0-1之间的随机标量】
% 初始化每个历史最优粒子
pBest = pos ; 
pbestfit = zeros(SearchAgents_no,1);
for i = 1:SearchAgents_no
pbestfit(i) = sum(pos(i,:).^2) ; 
end
%初始化全局历史最优粒子
[gBestfit,index] = min(pbestfit) ;
gBest = pos(index,:) ;
Convergence_curve = zeros(Max_iter,1);
 
for t=1:Max_iter
    for i=1:SearchAgents_no
        % 更新个体的位置和速度
        v(i,:) = w*v(i,:)+c1*rand*(pBest(i,:)-pos(i,:))+c2*rand*(gBest-pos(i,:)) ;
        pos(i,:) = pos(i,:)+v(i,:) ;
        % 边界处理
        v(i,:) = min(v(i,:), vmax);
        v(i,:) = max(v(i,:), -vmax);
        pos(i,:) =min(pos(i,:), ub);
        pos(i,:) =max(pos(i,:), lb);
        % 更新个体最优
        f1 = sum(pos(i,:).^2);
        if f1<pbestfit(i)    
           pBest(i,:) = pos(i,:) ;
           pbestfit(i) = f1;
        end
        % 更新全局最优
       if pbestfit(i) < gBestfit
            gBest = pBest(i,:) ;
            gBestfit = pbestfit(i) ;
       end
    end
    % 每代最优解对应的目标函数值
    Convergence_curve(t) = gBestfit; 
    disp(['Iteration = ' num2str(t)  ', Evaluations = ' num2str(gBestfit)]);
end
 
figure('unit','normalize','Position',[0.3,0.35,0.4,0.35],'color',[1 1 1],'toolbar','none')
subplot(1,2,1);
x = -5:0.1:5;y=x;
L=length(x);
f=zeros(L,L);
for i=1:L
    for j=1:L
       f(i,j) = x(i)^2+y(j)^2;
    end
end
surfc(x,y,f,'LineStyle','none');
xlabel('x_1');
ylabel('x_2');
zlabel('F')
title('Objective space')
 
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r','linewidth',1.5)
title('Convergence_curve')
xlabel('Iteration');
ylabel('Best score obtained so far');
 
axis tight
grid on
box on
legend('PSO')
display(['The best solution obtained by PSO is : ', num2str(gBest)]);
display(['The best optimal value of the objective funciton found by PSO is : ', num2str(gBestfit)]);
 
%% 
data=[4,4,5,5,5,6,6,6,6,7,7]
 [h,edges]=histcounts(data);
                    [~,idx]=max(h);
                    [~,idy]=min(h);
                    num_mode_max=h(idx);
                    num_mode_min=h(idy);
 disp(num_mode_max);
 disp(num_mode_min);

 %% 皮尔逊相关系数
 % 定义数据
temperature = [5, 7, 10, 15, 20, 25, 30, 28, 22, 16, 10, 6];
ice_cream_sales = [500, 600, 750, 850, 950, 1200, 1400, 1300, 1000, 900, 700, 550];

% 计算皮尔逊相关系数
correlation = corr(temperature', ice_cream_sales');

% 显示结果
disp(['气温和冰淇淋销量之间的皮尔逊相关系数为: ', num2str(correlation)]);

%% 多个变量的相关性分析
% 定义原始数据
data1 = randn(50, 1); % 第一组数据，50个观测值
data2 = randn(50, 1); % 第二组数据，50个观测值
data3 = randn(50, 1); % 第三组数据，50个观测值

% 取负值
data1_neg = -data1;
data2_neg = -data2;
data3_neg = -data3;

% 将负值数据合并为一个矩阵
data_neg = [data1,data2,data3,data1_neg, data2_neg, data3_neg];

% 计算相关矩阵
correlationMatrix_neg = corr(data_neg);

% 显示相关矩阵
disp('取负值后变量之间的相关矩阵为:');
disp(correlationMatrix_neg);

% 可视化相关矩阵
figure;
imagesc(correlationMatrix_neg);
colorbar;
title('Correlation Matrix (Negative Values)');
xlabel('Variables');
ylabel('Variables');

% 设置色彩映射
colormap('jet');

% 添加数值标签
numVariables = size(data_neg, 2);
for i = 1:numVariables
    for j = 1:numVariables
        text(j, i, num2str(correlationMatrix_neg(i, j), '%.2f'), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'black');
    end
end

% %% 3个数组的相关性分析
% % 示例数据
% A = [1, -2, 3, -4, 5];
% B = [-1, 2, -3, 4, -5];
% C = [1, -1, 1, -1, 1];
% 
% % 取负数组
% A_neg = -A;
% B_neg = -B;
% C_neg = -C;
% 
% % 处理异号元素的函数
% 1filter_arrays = @(arr1, arr2) arr1 .* arr2 > 0;
% 
% % 初始相关性矩阵
% correlation_matrix = zeros(3);
% 
% % 计算相关性
% arrays = {A_neg, B_neg, C_neg};
% for i = 1:3
%     for j = 1:3
%         if i ~= j
%             mask = filter_arrays(arrays{i}, arrays{j});
%             filtered_arr1 = arrays{i}(mask);
%             filtered_arr2 = arrays{j}(mask);
%             % 计算相关系数矩阵
%             corr_matrix = corrcoef(filtered_arr1, filtered_arr2);
%             % 获取相关系数矩阵中的第一个元素作为相关系数
%             correlation_matrix(i, j) = corr_matrix(1, 2);
%         else
%             correlation_matrix(i, j) = 1; % 自相关为1
%         end
%     end
% end
% 
% % 绘制热力图
% figure;
% heatmap({'A', 'B', 'C'}, {'A', 'B', 'C'}, correlation_matrix, 'Colormap', 'cool', 'ColorbarVisible', 'on', 'CellLabelColor','none');
% title('Correlation Heatmap');

% %% 示例
% % 示例数据
% A = [1, -2, 3, -4, 5];
% B = [-1, 2, -3, 4, -5];
% C = [1, -1, 1, -1, 1];
% 
% % 取负数组
% A_neg = -A;
% B_neg = -B;
% C_neg = -C;
% 
% % 处理异号元素的函数
% filter_arrays = @(arr1, arr2) arr1 .* arr2 > 0;
% 
% % 初始相关性矩阵
% correlation_matrix = zeros(3);
% 
% % 计算相关性
% arrays = {A_neg, B_neg, C_neg};
% for i = 1:3
%     for j = 1:3
%         if i ~= j
%             mask = filter_arrays(arrays{i}, arrays{j});
%             filtered_arr1 = arrays{i}(mask);
%             filtered_arr2 = arrays{j}(mask);
%             correlation_matrix(i, j) = corrcoef(filtered_arr1, filtered_arr2);
%         else
%             correlation_matrix(i, j) = 1; % 自相关为1
%         end
%     end
% end
% 
% % 绘制热力图
% figure;
% heatmap({'A', 'B', 'C'}, {'A', 'B', 'C'}, correlation_matrix, 'Colormap', 'cool', 'ColorbarVisible', 'on', 'CellLabelColor','none');
% title('Correlation Heatmap');

% %% 
% % 示例数据
% A = [1, -2, 3, -4, 5];
% B = [-1, 2, -3, 4, -5];
% C = [1, -1, 1, -1, 1];
% 
% % 原数组
% arrays_original = {A, B, C};
% % 取负数组,元胞
% arrays_negative = cellfun(@(x) -x, arrays_original, 'UniformOutput', false);
% 
% % 初始化相关性矩阵
% correlation_matrix = zeros(3);
% 
% % 计算原数组和负数组之间的相关性
% for i = 1:3
%     for j = 1:3
%         [filtered_arr1, filtered_arr2] = filter_arrays(arrays_original{i}, arrays_negative{j});
%         correlation_matrix(i, j) = calculate_correlation(filtered_arr1, filtered_arr2);
%     end
% end
% 
% % % 绘制热力图
% % figure;
% % colormap('cool'); % 设置颜色映射
% % heatmap({'-A', '-B', '-C'}, {'A', 'B', 'C'}, correlation_matrix, 'Colormap', 'cool', 'ColorbarVisible', 'on', 'CellLabelColor', 'black');
% % title('Correlation Heatmap');
% 
% figure;
% h = heatmap({'A', 'B', 'C'}, {'A', 'B', 'C'}, correlation_matrix, 'ColorbarVisible', 'on', 'CellLabelColor', 'red');
% colormap('parula'); % 设置颜色映射
% title('Correlation Heatmap');
% 
% % 处理异号元素的函数
% function [filtered_arr1, filtered_arr2] = filter_arrays(arr1, arr2)
%     mask = (arr1 .* arr2) > 0; % 保留同号元素
%     filtered_arr1 = arr1(mask);
%     filtered_arr2 = arr2(mask);
% end
% 
% % 计算相关性的函数
% function r = calculate_correlation(arr1, arr2)
%     if isempty(arr1) || isempty(arr2)
%         r = 0;
%     else
%         r = corrcoef(arr1, arr2);
%         r = r(1, 2); % 取相关系数矩阵中的相关系数
%     end
% end

% %% penalty function罚函数
% function [x_opt, f_opt] = penalty_method(f, g, x0, rho, epsilon, max_iter)
% % 惩罚函数法求解最优化问题
% % 输入参数：
% %   - f: 目标函数
% %   - g: 约束函数
% %   - x0: 初始点
% %   - rho: 惩罚参数
% %   - epsilon: 收敛准则
% %   - max_iter: 最大迭代次数
% % 输出结果：
% %   - x_opt: 最优解
% %   - f_opt: 最优解对应的目标函数值
%  
% x = x0;
% iter = 0;
% while iter < max_iter
%     % 构建带惩罚项的目标函数
%     p = @(x) f(x) + (rho/2) * sum(max(0, g(x)).^2);
%     
%     % 使用无约束优化算法求解带惩罚项的问题
%     options = optimoptions('fminunc','Algorithm','quasi-newton');
%     [x_opt, f_opt] = fminunc(p, x, options);
%     
%     % 判断是否达到收敛准则 
%     if norm(g(x_opt), 'inf') <= epsilon
%         break;
%     end
%     
%     % 更新惩罚参数
%     rho = rho * 10;
%     
%     x = x_opt;
%     iter = iter + 1;
% end

%%
% 示例数据
clusters = { [1, 2, 3], [4, 5] }; % 建筑编号
connections = { [0 1 1; 1 0 1; 1 1 0], [0 1; 1 0] }; % 连接矩阵

% 计算总建筑数量
totalBuildings = max(cellfun(@max, clusters));

% 初始化大矩阵
bigMatrix = zeros(totalBuildings);

% 遍历每个集群
for i = 1:length(clusters)
    cluster = clusters{i};
    connMatrix = connections{i};
    
    % 确定集群的建筑编号
    numBuildings = length(cluster);
    
    % 更新大矩阵中的连接关系
    for j = 1:numBuildings
        for k = 1:numBuildings
            bigMatrix(cluster(j), cluster(k)) = connMatrix(j, k);
        end
    end
end

% 确保矩阵对称
bigMatrix = (bigMatrix + bigMatrix')/2;

%%
% 示例数据
clusters = { [1, 2, 3], [4, 5]}; % 建筑编号
connections = { [0 1 1; 1 0 1; 1 1 0], [0 1; 1 0] }; % 连接矩阵

% 计算总建筑数量
totalBuildings = 0;
for i = 1:length(clusters)
    totalBuildings = max(totalBuildings, max(clusters{i}));
end

% 初始化大矩阵
bigMatrix = zeros(totalBuildings);

% 遍历每个集群
for i = 1:length(clusters)
    cluster = clusters{i};
    connMatrix = connections{i};
    
    % 更新大矩阵中的连接关系
    for j = 1:length(cluster)
        for k = 1:length(cluster)
            bigMatrix(cluster(j), cluster(k)) = connMatrix(j, k);
        end
    end
end

% 确保矩阵对称
bigMatrix = bigMatrix + bigMatrix';


%% 绘图
% 假设的建筑坐标数据 (20x2 矩阵，每行一个建筑的 (x, y) 坐标)
coords = [rand(20,1)*10, rand(20,1)*10]; % 示例数据

% 假设的0-1连接矩阵 (20x20 矩阵)
connectMatrix = randi([0, 1], 20, 20); % 示例数据，实际应替换为你的矩阵

% 绘制建筑坐标
figure;
hold on;
plot(coords(:,1), coords(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
text(coords(:,1), coords(:,2), num2str((1:20)'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% 绘制连接线
[n, m] = size(connectMatrix);
for i = 1:n
    for j = i+1:m
        if connectMatrix(i,j) == 1
            % 绘制建筑i和建筑j之间的连线
            plot([coords(i,1) coords(j,1)], [coords(i,2) coords(j,2)], 'b-');
        end
    end
end

% 设置图形属性
xlabel('X坐标');
ylabel('Y坐标');
title('建筑及其连接关系');
grid on;
axis equal;
hold off;

%%
% 假设的建筑坐标数据 (20x2 矩阵，每行一个建筑的 (x, y) 坐标)
coords = [rand(20,1)*10, rand(20,1)*10]; % 示例数据

% 假设的0-1连接矩阵 (20x20 矩阵)
connectMatrix = randi([0, 1], 20, 20); % 示例数据，实际应替换为你的矩阵

% 绘制建筑坐标
figure;
hold on;

% 设置颜色和线条样式
markerColor = [0.8, 0.2, 0.2]; % 红色
markerSize = 8; % 标记大小
lineColor = [0.2, 0.6, 1.0]; % 蓝色
lineWidth = 2; % 线宽

% 绘制建筑位置
plot(coords(:,1), coords(:,2), 'o', 'MarkerEdgeColor', markerColor, 'MarkerFaceColor', markerColor, 'MarkerSize', markerSize, 'LineWidth', 1.5);
text(coords(:,1), coords(:,2), num2str((1:20)'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% 绘制连接线
[n, m] = size(connectMatrix);
for i = 1:n
    for j = i+1:m
        if connectMatrix(i,j) == 1
            % 绘制建筑i和建筑j之间的连线
            plot([coords(i,1) coords(j,1)], [coords(i,2) coords(j,2)], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        end
    end
end

% 设置图形属性
xlabel('X坐标');
ylabel('Y坐标');
title('建筑及其连接关系');
grid on;
axis equal;
hold off;


