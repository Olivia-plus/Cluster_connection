% % 假设有50个建筑和24个时间点
% N = 50;
% T = 24;
% distance_threshold = 5; % 距离阈值，超过该值的建筑不考虑连接
% complementarity_weight = 10; % 互补性权重，增强互补性的重要性
% distance_weight = 1; % 距离权重，降低距离的重要性
% 
% % 生成显著的建筑净负荷特征数据，确保建筑互补性更明显
% rng(1); % 设置随机种子以便结果可重复
% net_load = [randn(N/2, T) + 2; randn(N/2, T) - 2]; % 每个建筑的日净负荷曲线，确保部分建筑负荷较高，部分较低
% coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀
% 
% % 计算相关性系数矩阵
% correlation_matrix = corr(net_load');
% 
% % 计算建筑之间的距离矩阵，并应用距离阈值
% distances = pdist2(coords, coords);
% distances(distances > distance_threshold) = inf; % 超过距离阈值的设为无穷大
% 
% % 计算互补度函数
% calculate_complementarity = @(A, B) sum(max(0, A) .* max(0, B)) + sum(min(0, A) .* min(0, B));
% 
% % 粒子群优化参数
% num_particles = 30;
% num_iterations = 100;
% w = 0.5; % 惯性权重
% c1 = 1.5; % 个人加速常数
% c2 = 1.5; % 社会加速常数
% 
% % 初始化粒子位置和速度，使用随机初始化
% positions = randi([1, ceil(N/2)], num_particles, N);
% velocities = zeros(num_particles, N);
% personal_best_positions = positions;
% global_best_position = positions(1, :);
% 
% % 初始化适应度
% fitness = inf(num_particles, 1);
% personal_best_fitness = fitness;
% global_best_fitness = inf;
% 
% 
% % PSO 迭代
% objective_values = zeros(num_iterations, 1); % 用于存储每次迭代的最优适应度值
% for iter = 1:num_iterations
%     for i = 1:num_particles
%         % 计算当前粒子的适应度
%         fitness(i) = fitness_function(positions(i, :), net_load, distances, correlation_matrix, calculate_complementarity, complementarity_weight, distance_weight);
%         
%         % 更新个人最优
%         if fitness(i) < personal_best_fitness(i)
%             personal_best_fitness(i) = fitness(i);
%             personal_best_positions(i, :) = positions(i, :);
%         end
%         
%         % 更新全局最优
%         if fitness(i) < global_best_fitness
%             global_best_fitness = fitness(i);
%             global_best_position = positions(i, :);
%         end
%     end
%     
%     % 更新粒子速度和位置，使用相关性系数矩阵指导粒子下一次迭代方向
%     for i = 1:num_particles
%         r1 = rand;
%         r2 = rand;
%         velocities(i, :) = w * velocities(i, :) + ...
%                            c1 * r1 * (personal_best_positions(i, :) - positions(i, :)) + ...
%                            c2 * r2 * (global_best_position - positions(i, :)) + ...
%                            correlation_matrix(i, :);
%         
%         positions(i, :) = positions(i, :) + velocities(i, :);
%         
%         % 确保位置在有效范围内
%         positions(i, positions(i, :) < 1) = 1;
%         positions(i, positions(i, :) > ceil(N/2)) = ceil(N/2);
%         positions(i, :) = round(positions(i, :));
%     end
%     
%     % 记录当前迭代的全局最优适应度值
%     objective_values(iter) = -global_best_fitness;
%     fprintf('Iteration %d, Global Best Fitness: %.2f\n', iter, global_best_fitness);
% end
% 
% % 最优集群划分结果
% best_labels = global_best_position;
% num_clusters = max(best_labels);
% 
% % 计算每个集群的最小生成树并可视化
% figure;
% hold on;
% colors = lines(num_clusters);
% total_complementarity = 0;
% total_cost = 0;
% for cluster = 1:num_clusters
%     cluster_buildings = find(best_labels == cluster);
%     scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
%     
%     % 计算集群内的最小生成树
%     cluster_distances = distances(cluster_buildings, cluster_buildings);
%     G = graph(cluster_distances, 'upper');
%     T = minspantree(G);
%     plot(G, 'XData', coords(cluster_buildings, 1), 'YData', coords(cluster_buildings, 2), 'EdgeColor', colors(cluster, :), 'LineWidth', 2);
%     
%     % 计算集群的总互补度和总互联成本
%     [r, c] = find(triu(T.Edges.EndNodes));
%     for k = 1:length(r)
%         total_complementarity = total_complementarity + correlation_matrix(cluster_buildings(r(k)), cluster_buildings(c(k)));
%         total_cost = total_cost + distances(cluster_buildings(r(k)), cluster_buildings(c(k)));
%     end
% end
% title(sprintf('最优聚类个数: %d, 总互补度: %.2f, 总互联成本: %.2f', num_clusters, total_complementarity, total_cost));
% xlabel('X 坐标');
% ylabel('Y 坐标');
% grid on;
% hold off;
% 
% % 绘制适应度值随迭代次数的变化
% figure;
% plot(1:num_iterations, objective_values, '-o');
% title('适应度值随迭代次数的变化');
% xlabel('迭代次数');
% ylabel('适应度值');
% grid on;
% 
% % 定义适应度函数
% function fit = fitness_function(position, net_load, distances, correlation_matrix, calculate_complementarity, complementarity_weight, distance_weight)
%     num_clusters = max(position);
%     total_complementarity = 0;
%     total_cost = 0;
%     
%     for cluster = 1:num_clusters
%         cluster_buildings = find(position == cluster);
%         if numel(cluster_buildings) < 2
%             fit = inf;
%             return;
%         end
%         
%         % 计算集群的互补度和互联成本
%         for i = 1:numel(cluster_buildings)
%             for j = i+1:numel(cluster_buildings)
%                 if distances(cluster_buildings(i), cluster_buildings(j)) < inf
%                     complementarity = calculate_complementarity(net_load(cluster_buildings(i), :), net_load(cluster_buildings(j), :));
%                     cost = distances(cluster_buildings(i), cluster_buildings(j));
%                     total_complementarity = total_complementarity + correlation_matrix(cluster_buildings(i), cluster_buildings(j));
%                     total_cost = total_cost + cost;
%                 end
%             end
%         end
%     end
%     
%     % 适应度值为负的总互补度减去总互联成本，使用权重调整
%     fit = -(complementarity_weight * total_complementarity - distance_weight * total_cost);
% end
% 


%%
% 假设有50个建筑和24个时间点
N = 50;
T = 24;
distance_threshold = 1; % 距离阈值，超过该值的建筑不考虑连接
complementarity_weight = 10; % 互补性权重，增强互补性的重要性
distance_weight = 1; % 距离权重，降低距离的重要性

% 生成显著的建筑净负荷特征数据，确保建筑互补性更明显
rng(1); % 设置随机种子以便结果可重复
net_load = [randn(N/2, T) + 2; randn(N/2, T) - 2]; % 每个建筑的日净负荷曲线，确保部分建筑负荷较高，部分较低
coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀

% 计算相关性系数矩阵
correlation_matrix = corr(net_load');

% 计算建筑之间的距离矩阵，并应用距离阈值
distances = pdist2(coords, coords);
distances(distances > distance_threshold) = inf; % 超过距离阈值的设为无穷大

% 计算互补度函数
calculate_complementarity = @(A, B) sum(max(0, A) .* max(0, B)) + sum(min(0, A) .* min(0, B));

% 粒子群优化参数
num_particles = 30;
num_iterations = 100;
w = 0.5; % 惯性权重
c1 = 1.5; % 个人加速常数
c2 = 1.5; % 社会加速常数

% 初始化粒子位置和速度，使用随机初始化
positions = randi([1, ceil(N/2)], num_particles, N);
velocities = zeros(num_particles, N);
personal_best_positions = positions;
global_best_position = positions(1, :);

% 初始化适应度
fitness = inf(num_particles, 1);
personal_best_fitness = fitness;
global_best_fitness = inf;

% PSO 迭代
objective_values = zeros(num_iterations, 1); % 用于存储每次迭代的最优适应度值
for iter = 1:num_iterations
    for i = 1:num_particles
        % 计算当前粒子的适应度
        fitness(i) = fitness_function(positions(i, :), net_load, distances, correlation_matrix, calculate_complementarity, complementarity_weight, distance_weight, distance_threshold);
        
        % 更新个人最优
        if fitness(i) < personal_best_fitness(i)
            personal_best_fitness(i) = fitness(i);
            personal_best_positions(i, :) = positions(i, :);
        end
        
        % 更新全局最优
        if fitness(i) < global_best_fitness
            global_best_fitness = fitness(i);
            global_best_position = positions(i, :);
        end
    end
    
    % 更新粒子速度和位置，使用相关性系数矩阵指导粒子下一次迭代方向
    for i = 1:num_particles
        r1 = rand;
        r2 = rand;
        velocities(i, :) = w * velocities(i, :) + ...
                           c1 * r1 * (personal_best_positions(i, :) - positions(i, :)) + ...
                           c2 * r2 * (global_best_position - positions(i, :)) + ...
                           correlation_matrix(i, :);
        
        positions(i, :) = positions(i, :) + velocities(i, :);
        
        % 确保位置在有效范围内
        positions(i, positions(i, :) < 1) = 1;
        positions(i, positions(i, :) > ceil(N/2)) = ceil(N/2);
        positions(i, :) = round(positions(i, :));
    end
    
    % 记录当前迭代的全局最优适应度值
    objective_values(iter) = -global_best_fitness;
    fprintf('Iteration %d, Global Best Fitness: %.2f\n', iter, global_best_fitness);
end

% 最优集群划分结果
best_labels = global_best_position;
num_clusters = max(best_labels);

% 计算每个集群的最小生成树并可视化
figure;
hold on;
colors = lines(num_clusters);% 给每个集群分配一种颜色
total_complementarity = 0;
total_cost = 0;
for cluster = 1:num_clusters
    cluster_buildings = find(best_labels == cluster);
    scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
    
    % 计算集群内的最小生成树
    cluster_distances = distances(cluster_buildings, cluster_buildings);
    G = graph(cluster_distances, 'upper');
    T = minspantree(G);
    plot(G, 'XData', coords(cluster_buildings, 1), 'YData', coords(cluster_buildings, 2), 'EdgeColor', colors(cluster, :), 'LineWidth', 2);
    
    % 计算集群的总互补度和总互联成本//1
    [r, c] = find(triu(T.Edges.EndNodes));
    for k = 1:length(r)
        total_complementarity = total_complementarity + correlation_matrix(cluster_buildings(r(k)), cluster_buildings(c(k)));
        total_cost = total_cost + distances(cluster_buildings(r(k)), cluster_buildings(c(k)));
    end
end
title(sprintf('最优聚类个数: %d, 总互补度: %.2f, 总互联成本: %.2f', num_clusters, total_complementarity, total_cost));
xlabel('X 坐标');
ylabel('Y 坐标');
grid on;
hold off;

% 绘制适应度值随迭代次数的变化
figure;
plot(1:num_iterations, objective_values, '-o');
title('适应度值随迭代次数的变化');
xlabel('迭代次数');
ylabel('适应度值');
grid on;


% 定义适应度函数
function fit = fitness_function(position, net_load, distances, correlation_matrix, calculate_complementarity, complementarity_weight, distance_weight, distance_threshold)
    num_clusters = max(position);
    total_complementarity = 0;
    total_cost = 0;
    modularity=0;
    
    for cluster = 1:num_clusters
        cluster_buildings = find(position == cluster);
        if numel(cluster_buildings) < 2
            fit = inf;
            return;
        end
        
        % 计算集群的互补度和互联成本
        for i = 1:numel(cluster_buildings)
            for j = i+1:numel(cluster_buildings)
                if distances(cluster_buildings(i), cluster_buildings(j)) < distance_threshold
                    complementarity = calculate_complementarity(net_load(cluster_buildings(i), :), net_load(cluster_buildings(j), :));
                    cost = distances(cluster_buildings(i), cluster_buildings(j));
                    total_complementarity = total_complementarity + correlation_matrix(cluster_buildings(i), cluster_buildings(j));
                    total_cost = total_cost + cost;
                else
                    fit = inf; % 如果有建筑距离超过阈值，则此解无效
                    return;
                end
            end
        end
    end
    
    % 适应度值为负的总互补度减去总互联成本，使用权重调整
    fit = -(complementarity_weight * total_complementarity - distance_weight * total_cost+modularity);
end
