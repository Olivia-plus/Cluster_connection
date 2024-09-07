% %% 1
% % 计算相似性矩阵
% alpha = 1.0;
% beta = 0.1;
% similarity_matrix = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         if i ~= j
%             complementarity = calculate_complementarity(net_load(i,:), net_load(j,:));
%             cost = distances(i, j);
%             similarity_matrix(i, j) = alpha * complementarity - beta * cost;
%         end
%     end
% end
% % 假设有N个建筑和T个时间点
% N = 20; % 建筑数量
% T = 48; % 时间点数量
% 
% % 生成示例数据
% net_load = randn(N, T); % 每个建筑的日净负荷曲线
% distances = rand(N, N); % 建筑之间的距离矩阵
% 
% % 使用谱聚类
% num_clusters = 3; % 设定集群数量
% L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
% [eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
% [~, idx] = sort(diag(eigValues)); % 按升序排序特征值
% U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）
% 
% % 进行K均值聚类
% labels = kmeans(U, num_clusters);
% 
% % 输出聚类结果
% for cluster = 1:num_clusters
%     disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
% end
% 
% 
% % 计算互补度函数
% function complementarity = calculate_complementarity(A, B)
%     complementarity = sum(max(0, A .* B));
% end

% %% 2
% % 假设有10个建筑和24个时间点
% N = 20;
% T = 48;
% 
% % 生成示例数据
% rng(1); % 设置随机种子以便结果可重复
% net_load = randn(N, T); % 每个建筑的日净负荷曲线
% distances = rand(N, N); % 建筑之间的距离矩阵
% distances = (distances + distances') / 2; % 确保矩阵对称
% 
% % 计算互补度函数
% calculate_complementarity = @(A, B) sum(max(0, A .* B));
% 
% % 计算相似性矩阵
% alpha = 1.0;
% beta = 0.1;
% similarity_matrix = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         if i ~= j
%             complementarity = calculate_complementarity(net_load(i,:), net_load(j,:));
%             cost = distances(i, j);
%             similarity_matrix(i, j) = alpha * complementarity - beta * cost;
%         end
%     end
% end
% 
% % 使用谱聚类
% num_clusters = 3;
% L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
% [eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
% [~, idx] = sort(diag(eigValues)); % 按升序排序特征值
% U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）
% 
% % 进行K均值聚类
% labels = kmeans(U, num_clusters);
% 
% % 绘制聚类结果
% figure;
% scatter3(U(:,1), U(:,2), U(:,3), 50, labels, 'filled');
% title('Spectral Clustering Result');
% xlabel('Eig1');
% ylabel('Eig2');
% zlabel('Eig3');
% grid on;
% 
% % 输出聚类结果
% for cluster = 1:num_clusters
%     disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
% end

% %% 3
% % 假设有10个建筑和24个时间点
% N = 10;
% T = 24;
% 
% % 生成示例数据，确保建筑均匀分布
% rng(1); % 设置随机种子以便结果可重复
% net_load = randn(N, T); % 每个建筑的日净负荷曲线
% coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀
% distances = squareform(pdist(coords)); % 计算建筑之间的距离矩阵
% 
% % 计算互补度函数
% calculate_complementarity = @(A, B) sum(max(0, A .* B));
% 
% % 计算相似性矩阵
% alpha = 1.0;
% beta = 0.1;
% similarity_matrix = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         if i ~= j
%             complementarity = calculate_complementarity(net_load(i,:), net_load(j,:));
%             cost = distances(i, j);
%             similarity_matrix(i, j) = alpha * complementarity - beta * cost;
%         end
%     end
% end
% 
% % 使用谱聚类
% num_clusters = 3;
% L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
% [eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
% [~, idx] = sort(diag(eigValues)); % 按升序排序特征值
% U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）
% 
% % 进行K均值聚类，确保均衡的集群划分
% options = statset('Display','final');
% labels = kmeans(U, num_clusters, 'Options', options, 'MaxIter', 1000);
% 
% % 可视化结果
% figure;
% hold on;
% colors = lines(num_clusters);
% for cluster = 1:num_clusters
%     cluster_buildings = find(labels == cluster);
%     scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
%     for i = 1:length(cluster_buildings)
%         for j = i+1:length(cluster_buildings)
%             if similarity_matrix(cluster_buildings(i), cluster_buildings(j)) > 0
%                 plot([coords(cluster_buildings(i), 1), coords(cluster_buildings(j), 1)], ...
%                      [coords(cluster_buildings(i), 2), coords(cluster_buildings(j), 2)], ...
%                      'Color', colors(cluster, :));
%             end
%         end
%     end
% end
% title('建筑的集群划分与互联情况');
% xlabel('X 坐标');
% ylabel('Y 坐标');
% grid on;
% hold off;
% 
% % 输出聚类结果
% for cluster = 1:num_clusters
%     disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
% end
% 
% %% 4
% % 假设有10个建筑和24个时间点
% N = 12;
% T = 24;
% 
% % 生成示例数据，确保建筑均匀分布
% rng(1); % 设置随机种子以便结果可重复
% net_load = randn(N, T); % 每个建筑的日净负荷曲线
% coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀
% distances = squareform(pdist(coords)); % 计算建筑之间的距离矩阵
% 
% % 计算互补度函数
% calculate_complementarity = @(A, B) sum(max(0, A .* B));
% 
% % 计算相似性矩阵
% alpha = 1.0;
% beta = 0.1;
% similarity_matrix = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         if i ~= j
%             complementarity = calculate_complementarity(net_load(i,:), net_load(j,:));
%             cost = distances(i, j);
%             similarity_matrix(i, j) = alpha * complementarity - beta * cost;
%         end
%     end
% end
% 
% % 使用谱聚类
% num_clusters = 5;
% L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
% [eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
% [~, idx] = sort(diag(eigValues)); % 按升序排序特征值
% U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）
% 
% % 进行K均值聚类，确保均衡的集群划分
% options = statset('Display','final');
% labels = kmeans(U, num_clusters, 'Options', options, 'MaxIter', 1000);
% 
% % 可视化结果
% figure;
% hold on;
% colors = lines(num_clusters);
% for cluster = 1:num_clusters
%     cluster_buildings = find(labels == cluster);
%     scatter(coords(cluster_buildings, 1), coords(cluster_buildings, 2), 100, 'filled', 'MarkerFaceColor', colors(cluster, :));
%     for i = 1:length(cluster_buildings)
%         for j = i+1:length(cluster_buildings)
%             if similarity_matrix(cluster_buildings(i), cluster_buildings(j)) > 0
%                 plot([coords(cluster_buildings(i), 1), coords(cluster_buildings(j), 1)], ...
%                      [coords(cluster_buildings(i), 2), coords(cluster_buildings(j), 2)], ...
%                      'Color', colors(cluster, :));
%             end
%         end
%     end
% end
% title('建筑的集群划分与互联情况');
% xlabel('X 坐标');
% ylabel('Y 坐标');
% grid on;
% hold off;
% 
% % 输出聚类结果
% for cluster = 1:num_clusters
%     disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
% end



% 
%% 5
% 假设有10个建筑和24个时间点
N = 10;
T = 24;

% 生成示例数据，确保建筑均匀分布
rng(1); % 设置随机种子以便结果可重复
net_load = randn(N, T); % 每个建筑的日净负荷曲线
coords = rand(N, 2) * 10; % 建筑的二维坐标，确保分布均匀
distances = squareform(pdist(coords)); % 计算建筑之间的距离矩阵

% 计算互补度函数
calculate_complementarity = @(A, B) sum(max(0, A .* B));

% 计算相似性矩阵
alpha = 1.0;
beta = 0.1;
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

% 使用谱聚类
num_clusters = 3;
L = diag(sum(similarity_matrix, 2)) - similarity_matrix; % 拉普拉斯矩阵
[eigVectors, eigValues] = eig(L); % 计算特征向量和特征值
[~, idx] = sort(diag(eigValues)); % 按升序排序特征值
U = eigVectors(:, idx(2:num_clusters+1)); % 选取最小的特征向量（不包括第一个）

% 进行K均值聚类，确保均衡的集群划分
options = statset('Display','final');
labels = kmeans(U, num_clusters, 'Options', options, 'MaxIter', 1000);

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
title(sprintf('总互补度: %.2f, 总互联成本: %.2f', total_complementarity, total_cost));
xlabel('X 坐标');
ylabel('Y 坐标');
grid on;
hold off;

% 输出聚类结果
for cluster = 1:num_clusters
    disp(['Cluster ' num2str(cluster) ': ' num2str(find(labels == cluster)')]);
end
