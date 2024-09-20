%**************************************************************************
% 文件名: C:\Users\WuAoli\Desktop\毕设\集群规划算法代码\Cluster\main.m
% 版本: v1.0
% 作者: Wuaoli-123
% 联系方式: 2713337051@qq.com
% 日期: 2024-09-08
% 描述:适应度函数由成本和收益两个指标构成。收益由集群整体构成的新增光伏消纳量决定，成本由线路铺设投资成本决定。去掉模块度，同时在集群划分的过程中，将储能和可平移负荷的综合优化融合到目标函数最优化的求解中去。
% 输入:  
% 输出:  
%**************************************************************************
clear;
clc;
format short;% 默认精度
tic
%% 设定建筑的数量，地理坐标，以及建筑类型的分类
% 设定建筑数量 20~100不等
num_buildings = 20;
% 随机生成建筑横、纵坐标；建筑类型；负荷和光伏曲线分配
% 建筑类型确定，负荷和光伏的曲线也就确定了【待优化】
[x, y, type, load_curve, pv_curve,flexible_load_main,storage_capacity_main] = GenerateBuildingInfo(num_buildings);
% 绘制建筑位置分布图
PlotBuildingLocations(x, y,type);
% 显示每个建筑的位置和类型, '行列转换符,列转行
building_info= table((1:num_buildings)', x', y', type', 'VariableNames', {'建筑编号', '横坐标', '纵坐标', '建筑类型'});
disp('随机生成的建筑位置和类型：');
disp(building_info);
% 绘制每个建筑的负荷曲线和光伏曲线以及净负荷图【最后再取消注释，这一步其实也不是很有必要，但是也写着】
% 循环绘制每个建筑的曲线图
% figure;
% for i = 1:num_buildings
%     % 绘制第i个建筑的曲线图
%     PlotBuildingCurves(load_curve{i}, pv_curve{i});
%     hold on
% end
% % 添加标题
% title('20个建筑曲线图');
% %% 绘制负荷变化三维图【不是很重要】
% %  净负荷元胞
% hours_per_day = 48;
% % 创建一个空的元胞数组来存储净负荷曲线
% net_load_curve = cell(1, num_buildings);
% % 净负荷       
% for i = 1:num_buildings
%     net_load_curve{i} = load_curve{i} - pv_curve{i};
% end
% % 创建时间轴标签，感觉没起到作用！！！【为啥捏？】
% time_labels = cell(1, hours_per_day);
% for i = 1:hours_per_day
%     hour = floor((i-1)/2);
%     minute = rem(i-1, 2) * 30;
%     time_labels{i} = sprintf('%02d:%02d', hour, minute);
% end
% % 创建小时和建筑索引网格
% [hour, building] = meshgrid(1:hours_per_day, 1:num_buildings);
% 
% % 将每个建筑的负荷数据放在网格中的对应位置
% load_grid = zeros(num_buildings, hours_per_day);
% for i = 1:num_buildings
%     load_grid(i, :) = net_load_curve{i};
% end
% 
% % 绘制三维图
% figure;
% surf(building, hour, load_grid);
% xlabel('建筑编号');
% ylabel('时间');
% % zlabel('负荷/MW');
% % zlabel('光伏/MW');
% zlabel('净负荷/MW');
% % title('所有建筑节点的日负荷数据');
% % title('所有建筑节点的日光伏出力数据');
% title('所有建筑节点的日净负荷数据');
% colorbar; % 添加颜色条
% grid on

%% 根据提供的数据编写粒子群算法的集群划分代码
%% 基本参数设置
max_iter = 50; % 最大迭代次数
pop_size = 100; % 种群规模
dim=num_buildings; % 粒子维度
% ub=num_buildings; %集群划分个数上限
% lb=1; %集群划分个数下限
vmax = 4; % 最大飞行速度
max_num_cluster =ceil(num_buildings/2); % 最大集群划分数量为建筑的总个数/2，ceil向上取整
w = 0.8; % 惯性权重
c1 = 1.5; % 学习因子 1
c2 = 1.5; % 学习因子 2        
electricity_price=0.45;% 建筑交易收益电价0.45元/度，恒定不变
dc_cost_p=208;% 线路铺设成本208元/米【TODO:折旧到每一天】

net_load{num_buildings}=0;
    for i = 1:num_buildings
    net_load{i} = load_curve{i} - pv_curve{i};% 净负荷曲线
    end
% 初始化种群
particles = zeros(pop_size, dim); % 初始化种群的位置 每行代表一个粒子集群划分情况，每列代表一个建筑，值表示所属集群编号 100行n列
velocity=zeros(pop_size, dim);
% 初始化每个历史最优粒子
pbest_fitness = Inf(pop_size, 1); % 个体历史最优适应度值  总成本最低
trade_power_1 = zeros(pop_size, 1);% 交易的电量
% bigMatrix_1=zeros(pop_size, 1);
% 初始化粒子位置，可行解就是建筑随机分到种群的过程，因为n个建筑随机分到多个种群容器中会有很多种情况，遍历起来会很麻烦。
% 粒子信息初始化
for i = 1:pop_size
    particles(i, :) = randi(max_num_cluster, 1, num_buildings);
    velocity(i,:) = -vmax+2*vmax*randi(vmax,1,num_buildings); % 初始化速度
%             % 检查这个随机解中是否每个集群包含至少一个建筑【感觉没必要】
%             while ~any(ismember(1:num_buildings, particles(i, :)))
%             % 如果没有建筑，重新生成一个随机解
%             particles(i, :) = randi(p, 1, num_buildings);
%             end
    [pbest_fitness(i,1),trade_power_1(i,1),]=calculate_fitness(particles(i, :),net_load, electricity_price, dc_cost_p,x,y,num_buildings,flexible_load_main,storage_capacity_main);
end
pbest = particles; % 所有粒子个体最优位置  
% 初始化全局历史最优粒子 
[gbest_fitness,index]=min(pbest_fitness);
gbest=particles(index,:);
trade=0;
Convergence_curve=zeros(max_iter,1); % 收敛曲线
trade_curve=zeros(max_iter,1); % 交易曲线
best_connectMatrix=zeros(num_buildings,num_buildings);% 最佳连接矩阵
% for p=3:num_buildings % 这n个建筑可以划分为1~n个种群，只是用3来做测试
        % 迭代优化
        fitness_valuse_personal=zeros(pop_size,1);
        trade_power=zeros(pop_size,1);
        bigMatrix=zeros(num_buildings,num_buildings);
        % 迭代开始
        for iter = 1:max_iter 
            % 对所有的粒子遍历
            for j = 1:pop_size
                    % 更新个体的位置和速度
                    [particles(j,:),velocity(j,:)] = update_particle_position(particles(j, :), pbest(j, :), gbest, w, c1, c2,velocity(j,:),vmax,max_num_cluster); 
                    % 求最大集群中建筑的个数
                    [h,edges]=histcounts(particles(j, :));
                    [~,idx]=max(h);
                    [~,idy]=min(h);
                    num_mode_max=h(idx);
                    num_mode_min=h(idy);
                    % 做一个小的判断，只有满足划分要求的粒子才能进行适应度的计算
                   if num_mode_max < 8 && num_mode_min > 1
                        % 计算当前粒子的适应度值
                        [fitness_valuse_personal(j,1),trade_power(j,1),bigMatrix]= calculate_fitness(particles(j, :), net_load, electricity_price, dc_cost_p,x,y,num_buildings,flexible_load_main,storage_capacity_main); 
                        % 更新个体最优
                        if fitness_valuse_personal(j,1) < pbest_fitness(j,1)
                            pbest(j,:) = particles(j, :);
                            pbest_fitness(j,1) = fitness_valuse_personal(j,1);
                        end
                        % 更新全局最优
                        if pbest_fitness(j,1) < gbest_fitness% 粒子和全局最优解对比
                            gbest = pbest(j,:);
                            gbest_fitness = pbest_fitness(j,1);
                            trade=trade_power(j,1);
                            best_connectMatrix=bigMatrix;
                        end
                   else
                       fitness_valuse_personal(j,1)=inf;
                   end
            end
            % 每代最优解对应的目标函数值
            Convergence_curve(iter)=gbest_fitness;
            trade_curve(iter)=trade;
            disp(['Iteration = ' num2str(iter) ', Evaluations = ' num2str(gbest)]);
        end
        disp(['集群最优适应度为 = ' num2str(-Convergence_curve(max_iter)) '元， '' 集群光伏总消纳量为 = ' num2str(trade_curve(max_iter)/2) 'kWh']);
%         %% 显示最终结果
%         disp('最优集群划分方案：');
%         disp(gbest);
%         disp('全局最优适应度值：');
%         disp(gbest_fitness);
       
        figure 
        plot(1:max_iter, Convergence_curve);
        xlabel('迭代次数');
        ylabel('适应度函数值');
        title('适应度函数随迭代次数的变化');
         
%         figure
%         plot(1:max_iter, fitness_valuse_personal);
%         xlabel('迭代次数');
%         ylabel('适应度函数值');
%         title('粒子适应度函数随迭代次数的变化');
% end

% %% 绘制相关性热力图
% % 示例数据
% % A = [1, -2, 3, -4, 5];
% % B = [-1, 2, -3, 4, -5];
% % C = [1, -1, 1, -1, 1];
% 
% % 原数组
% arrays_original = net_load_curve;
% % 取负数组,元胞
% arrays_negative = cellfun(@(x) -x, arrays_original, 'UniformOutput', false);
% 
% % 初始化相关性矩阵
% correlation_matrix = zeros(num_buildings);
% 
% % 计算原数组和负数组之间的相关性
% for i = 1:num_buildings
%     for j = 1:num_buildings
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
% % h = heatmap({'1', '2', '3','1', '2', '3','1', '2', '3','1', '2', '3','1', '2', '3','1', '2', '3','1', '2'}, {'A', 'B', 'C','A', 'B', 'C','A', 'B', 'C','A', 'B', 'C','A', 'B', 'C','A', 'B', 'C','A', 'B'}, correlation_matrix, 'ColorbarVisible', 'on', 'CellLabelColor', 'red');
% h = heatmap( correlation_matrix, 'ColorbarVisible', 'on', 'CellLabelColor', 'black');
% % colormap('parula'); % 设置颜色映射
% colormap('default'); % 设置颜色映射
% title('Correlation Heatmap');
% 
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

% % 最优集群划分结果
% best_labels = gbest;
% num_clusters = max(best_labels);
% % 集群连接可视化
% figure;
% hold on;
% colors = lines(num_clusters);% 给每个集群分配一种颜色
% for cluster = 1:num_clusters
%     cluster_buildings = find(best_labels == cluster);
%     scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
%     
%     % 计算集群内的最小生成树
%     cluster_distances = distances(cluster_buildings, cluster_buildings);
%     G = graph(cluster_distances, 'upper');
%     T = minspantree(G);
%     plot(G, 'XData', coords(cluster_buildings, 1), 'YData', coords(cluster_buildings, 2), 'EdgeColor', colors(cluster, :), 'LineWidth', 2);
% end
% title(sprintf('最优聚类个数: %d', num_clusters));
% xlabel('X 坐标');
% ylabel('Y 坐标');
% grid on;
% hold off;

% 假设的建筑坐标数据 (每行一个建筑的 (x, y) 坐标)
coords = [x',y']; 

% 假设的0-1连接矩阵 (20x20 矩阵)
connectMatrix = best_connectMatrix; % 示例数据，实际应替换为你的矩阵
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
text(coords(:,1), coords(:,2), num2str((1:num_buildings)'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

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

