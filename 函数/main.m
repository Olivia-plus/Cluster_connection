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
% disp('随机生成的建筑位置和类型：');
% disp(building_info);
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
%                                            
% grid on


% kmeans算法，将100个粒子按照不同的种群的数量（随机给）进行聚类。
% 然后每个粒子按照给定的速度和方向去扩张和缩减自己所属种群的大小，但是仍旧保证集群之间不存在交叉重叠的情况。


%% 根据提供的数据编写粒子群算法的集群划分代码
% 基本参数设置
max_iter = 5; % 最大迭代次数
pop_size = 4; % 种群规模
dim=num_buildings; % 粒子维度
% max_num_cluster =ceil(num_buildings/4); % 最大集群划分数量为建筑的总个数/2，ceil向上取整
numClusters = 5;         % 集群数量

w = 0.5; % 惯性权重
c1 = 1.5; % 学习因子 1
c2 = 1.5; % 学习因子 2        
velocityLimit=100;% 粒子速度限制
coord=[x',y'];% 坐标
electricity_price=0.25;% 建筑交易收益电价0.25元/度，恒定不变

net_load{num_buildings}=0;
    for i = 1:num_buildings
    net_load{i} = load_curve{i} - pv_curve{i};% 净负荷曲线
    end

% 显示初始点分布
figure;
scatter(coord(:,1), coord(:,2), 'filled');
title('建筑分布');
xlabel('X');
ylabel('Y');
hold on;

% 初始化粒子群
particles = cell(pop_size, 1);  % 每个粒子表示不同的质心集
velocities = cell(pop_size, 1); % 每个粒子的速度

for i = 1:pop_size
    particles{i} = coord(randperm(num_buildings, numClusters), :);  % 随机选择质心
    velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
end

% 初始化每个历史最优粒子
pbest_fitness = Inf(pop_size, 1); 
pbest = particles; 
gbest_fitness=Inf;
gbest=particles{1};

% 记录适应度变化
fitnessHistory = zeros(max_iter, 1);  % 适应度变化记录

trade=0;
Convergence_curve=zeros(max_iter,1); % 收敛曲线
trade_curve=zeros(max_iter,1); % 交易曲线
best_connectMatrix=zeros(num_buildings,num_buildings);% 最佳连接矩阵

% for p=3:num_buildings % 这n个建筑可以划分为1~n个种群，只是用3来做测试
        fitness_valuse_personal=zeros(pop_size,1);
        trade_power=zeros(pop_size,1);
        bigMatrix=cell(1,pop_size);
       
        % 迭代开始
        for iter = 1:max_iter 
            % 对所有的粒子遍历
            for j = 1:pop_size
                % 获取当前粒子的质心
                C = particles{j};
                % 计算每个点到质心的距离
                distances = pdist2(coord, C);
                [~, clusterIdx] = min(distances, [], 2);  % 将每个点分配到最近的质心
                % 检查集群是否满足要求：没有孤立节点，集群不重叠
                valid = checkClusterValidity(coord, clusterIdx, numClusters);
                % 计算适应度
                if valid
                    [fitness_valuse_personal(j),trade_power(j),bigMatrix{j}]= calculate_fitness(clusterIdx,load_curve,pv_curve,electricity_price,x,y,num_buildings,flexible_load_main,storage_capacity_main); % 【将net_load替换成了load_curve,pv_curve,便于计算柔性负荷最优调度】
                else
                    fitness_valuse_personal(j) = inf;  % 不合法的集群划分给予惩罚
                end

                % 更新个体最优
                if fitness_valuse_personal(j) < pbest_fitness(j)
                    pbest{j}= particles{j};
                    pbest_fitness(j) = fitness_valuse_personal(j);
                end

                % 更新全局最优
                if pbest_fitness(j) < gbest_fitness% 粒子和全局最优解对比
                    gbest = particles{j};
                    gbest_fitness = fitness_valuse_personal(j);
                    trade=trade_power(j);
                    best_connectMatrix=bigMatrix{j};
                end
            end

            % 记录当前迭代的适应度
            fitnessHistory(iter) = gbest_fitness;
        
            % 更新质心位置
            for p = 1:pop_size
                velocities{p} = w * velocities{p} ...
                    + c1 * rand * (pbest{p} - particles{p}) ...
                    + c2 * rand * (gbest - particles{p});  % 粒子速度更新
                particles{p} = particles{p} + velocities{p};  % 更新粒子质心位置
            end

            % 可视化第一个粒子的集群划分
            figure(1);
            subplot(1, 2, 1);
            visualizeClusters(coord, particles{1}, clusterIdx, numClusters);
            title(['第一个粒子的集群划分 (迭代: ', num2str(iter), ')']);

            % 可视化当前最优集群划分（每隔10次迭代）
            subplot(1, 2, 2);
            visualizeClusters(coord, gbest, clusterIdx, numClusters);
            title(['全局最优集群划分 (迭代: ', num2str(iter), ')']);

            pause(0.5);

            % 检查是否达到收敛条件
            if gbest_fitness < 1e-5
                break;
            end
        end

        % 显示最终集群划分结果
        figure(2);
        visualizeClusters(coord, gbest, clusterIdx, numClusters);
        title('最终最优集群划分');
        
        % 绘制适应度变化图
        figure(3);
        plot(fitnessHistory(1:iter), 'LineWidth', 2);
        title('适应度变化图');
        xlabel('迭代次数');
        ylabel('适应度');
        grid on;

        disp(['集群最优适应度为 = ' num2str(-gbest) '元， '' 集群光伏总消纳量为 = ' num2str(trade) 'kWh']);

%% 绘制建筑互联图
connectMatrix = best_connectMatrix; % 示例数据，实际应替换为你的矩阵
num_buildings = size(coords, 1); % 建筑数量

% 获取连通组件（建筑群体）
G = graph(connectMatrix); % 将连接矩阵转为图
[bin, binsizes] = conncomp(G); % bin表示每个节点所属的连通分量，binsizes表示每个分量的大小

% 获取颜色
unique_bins = unique(bin);
num_clusters = length(unique_bins); % 连通子图数量
colors = lines(num_clusters); % 使用不同颜色表示不同的连通子图

% 设置线条样式
lineColor_other = [0.5, 0.5, 0.5]; % 灰色，用于不同连通分量之间的线
lineWidth = 2; % 线宽
markerSize = 8; % 标记大小

% 绘制建筑位置
figure;
hold on;

% 绘制所有建筑的节点
plot(coords(:,1), coords(:,2), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', markerSize, 'LineWidth', 1.5);
text(coords(:,1), coords(:,2), num2str((1:num_buildings)'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% 绘制连通的树枝和节点圈
for k = 1:num_clusters
    % 获取属于当前连通分量的建筑索引
    cluster_nodes = find(bin == unique_bins(k));
    
    % 如果连通分量内的建筑数量大于2，画圈圈出这些节点
    if length(cluster_nodes) > 2
        cluster_coords = coords(cluster_nodes, :);
        hull = convhull(cluster_coords(:,1), cluster_coords(:,2)); % 获取凸包
        fill(cluster_coords(hull,1), cluster_coords(hull,2), colors(k,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 轻微填充
    elseif length(cluster_nodes) == 2
        % 当群体只有两个节点时，用矩形圈住
        cluster_coords = coords(cluster_nodes, :);
        rectangle('Position', [min(cluster_coords(:,1)) min(cluster_coords(:,2)) ...
            abs(diff(cluster_coords(:,1))) abs(diff(cluster_coords(:,2)))], ...
            'EdgeColor', colors(k,:), 'LineWidth', 1.5, 'LineStyle', '--');
    end
    
    % 为当前群体绘制线（包括仅两个节点的群体）
    for i = 1:length(cluster_nodes)
        for j = i+1:length(cluster_nodes)
            node1 = cluster_nodes(i);
            node2 = cluster_nodes(j);
            if connectMatrix(node1, node2) == 1
                % 绘制建筑 node1 和 node2 之间的连线，使用特定颜色
                plot([coords(node1,1) coords(node2,1)], [coords(node1,2) coords(node2,2)], '-', 'Color', colors(k,:), 'LineWidth', lineWidth);
            end
        end
    end
end

% 绘制其他未连通的建筑之间的线
[n, m] = size(connectMatrix);
for i = 1:n
    for j = i+1:m
        if connectMatrix(i,j) == 1 && bin(i) ~= bin(j)
            % 绘制不同连通分量间的连线，使用灰色
            plot([coords(i,1) coords(j,1)], [coords(i,2) coords(j,2)], '--', 'Color', lineColor_other, 'LineWidth', lineWidth);
        end
    end
end

hold off;


% 设置图形属性
xlabel('X坐标');
ylabel('Y坐标');
title('建筑及其连接关系');
grid on;
axis equal;
hold off;

%% 辅助函数：检查集群有效性
function valid = checkClusterValidity(X, clusterIdx, numClusters)
    for i = 1:numClusters
        if sum(clusterIdx == i) < 2
            valid = false;
            return;
        end
    end
    
    clusterCenters = zeros(numClusters, 2);
    for i = 1:numClusters
        clusterCenters(i, :) = mean(X(clusterIdx == i, :), 1);
    end
    
    if min(pdist(clusterCenters)) < 50
        valid = false;
        return;
    end
    valid = true;
end

%% 辅助函数：可视化集群划分
function visualizeClusters(X, C, clusterIdx, numClusters)
    colors = lines(numClusters);  % 使用不同的颜色标识不同集群
    hold on;
    
    for i = 1:numClusters
        clusterPoints = X(clusterIdx == i, :);  % 选出属于第 i 个集群的点
        scatter(clusterPoints(:, 1), clusterPoints(:, 2), 50, colors(i, :), 'filled');
        
        % 画出集群的边界（凸包）
        if size(clusterPoints, 1) > 2
            k = convhull(clusterPoints(:, 1), clusterPoints(:, 2));
            plot(clusterPoints(k, 1), clusterPoints(k, 2), 'Color', colors(i, :), 'LineWidth', 2);
        end
        
        scatter(C(i, 1), C(i, 2), 100, 'x', 'LineWidth', 2, 'MarkerEdgeColor', colors(i, :));  % 质心
    end
    
    title('集群划分结果');
    xlabel('X');
    ylabel('Y');
    hold off;
end
