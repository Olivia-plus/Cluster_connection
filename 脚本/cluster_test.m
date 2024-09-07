% 假设有10个建筑和24个时间点
N = 16;
T = 48;

% 生成显著的建筑净负荷特征数据，确保建筑互补性更明显
rng(1); % 设置随机种子以便结果可重复
net_load = [randn(N/2, T) + 10; randn(N/2, T) - 10]; % 每个建筑的日净负荷曲线，确保部分建筑负荷较高，部分较低
coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀
distances = squareform(pdist(coords)); % 计算建筑之间的距离矩阵

% 计算互补度函数
calculate_complementarity = @(A, B) sum(max(0, A) .* max(0, B)) + sum(min(0, A) .* min(0, B));

% 计算相似性矩阵
alpha = 10.0; % 加大互补性的权重
beta = 1.0;
similarity_matrix = zeros(N, N);
for i = 1:N
    for j = 1:N
        if i ~= j
            complementarity = calculate_complementarity(net_load(i,:), net_load(j,:));
            cost = distances(i, j);
            similarity_matrix(i, j) = alpha * complementarity - beta * cost;
        end
    end
end

% 寻找最优的聚类个数
best_silhouette = -1;
best_num_clusters = 1;
best_labels = [];
objective_values = [];

for num_clusters = 2:N-1
    L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
    [eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
    [~, idx] = sort(diag(eigValues)); % 按升序排序特征值
    U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）
    
    % 进行K均值聚类
    labels = kmeans(U, num_clusters, 'MaxIter', 1000);
    
    % 检查是否存在孤立的集群
    valid_clustering = true;
    for cluster = 1:num_clusters
        if sum(labels == cluster) < 2
            valid_clustering = false;
            break;
        end
    end
    
    % 如果存在孤立集群，则跳过该聚类个数
    if ~valid_clustering
        continue;
    end
    
    % 计算轮廓系数
    silhouette_values = silhouette(coords, labels);
    mean_silhouette = mean(silhouette_values);
    
    % 计算集群的目标函数值
    total_complementarity = 0;
    total_cost = 0;
    for cluster = 1:num_clusters
        cluster_buildings = find(labels == cluster);
        
        % 计算集群内的最小生成树
        cluster_distances = distances(cluster_buildings, cluster_buildings);
        G = graph(cluster_distances, 'upper');
        T = minspantree(G);
        
        % 计算集群的总互补度和总互联成本
        [r, c] = find(triu(T.Edges.EndNodes));
        for k = 1:length(r)
            total_complementarity = total_complementarity + similarity_matrix(cluster_buildings(r(k)), cluster_buildings(c(k)));
            total_cost = total_cost + distances(cluster_buildings(r(k)), cluster_buildings(c(k)));
        end
    end
    objective_value = total_complementarity - total_cost;
    objective_values = [objective_values; num_clusters, objective_value];
    
    % 如果当前轮廓系数更高，则更新最优聚类个数
    if mean_silhouette > best_silhouette
        best_silhouette = mean_silhouette;
        best_num_clusters = num_clusters;
        best_labels = labels;
    end
end

% 确保 objective_values 非空
if isempty(objective_values)
    error('没有找到有效的聚类结果，请检查输入数据或参数设置。');
end

% 使用最优的聚类个数和标签
labels = best_labels;
num_clusters = best_num_clusters;

% 计算每个集群的最小生成树并可视化
figure;
hold on;
colors = lines(num_clusters);
total_complementarity = 0;
total_cost = 0;
for cluster = 1:num_clusters
    cluster_buildings = find(labels == cluster);
    scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
    
    % 计算集群内的最小生成树
    cluster_distances = distances(cluster_buildings, cluster_buildings);
    G = graph(cluster_distances, 'upper');
    T = minspantree(G);
    
    % 绘制最小生成树
    [r, c] = find(triu(T.Edges.EndNodes));
    for k = 1:length(r)
        plot([coords(cluster_buildings(r(k)), 1), coords(cluster_buildings(c(k)), 1)], ...
             [coords(cluster_buildings(r(k)), 2), coords(cluster_buildings(c(k)), 2)], ...
             'Color', colors(cluster, :));
        total_complementarity = total_complementarity + similarity_matrix(cluster_buildings(r(k)), cluster_buildings(c(k)));
        total_cost = total_cost + distances(cluster_buildings(r(k)), cluster_buildings(c(k)));
    end
end
title(sprintf('最优集群划分 (集群数: %d) - 总互补度: %.2f, 总互联成本: %.2f', num_clusters, total_complementarity, total_cost));
xlabel('X 坐标');
ylabel('Y 坐标');
grid on;
hold off;

% 输出聚类结果
for cluster = 1:num_clusters
    disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
end

% 可视化目标函数值随着迭代次数的变化
figure;
plot(objective_values(:,1), objective_values(:,2), '-o');
xlabel('集群个数');
ylabel('目标函数值 (总互补度 - 总互联成本)');
title('目标函数值随着集群个数的变化');
grid on;


