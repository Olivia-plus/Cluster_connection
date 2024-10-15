% % 示例数据
% data = [
%     1, 2;
%     1.5, 1.8;
%     5, 8;
%     8, 8;
%     1, 0.6;
%     9, 11
% ];
%  
% % 设置 DBSCAN 参数
% epsilon = 2;
% minPts = 2;
%  
% % 运行 DBSCAN
% [clusters, clusterCount] = dbscan(data, epsilon, minPts);
%  
% % 打印聚类结果
% disp('聚类结果：');
% disp(clusters);
% disp('聚类数目：');
% disp(clusterCount);
% % DBSCAN算法
% function [clusters, clusterCount] = dbscan(data, epsilon, minPts)
%     dataSize = size(data, 1);
%     visited = false(dataSize, 1);
%     clusters = zeros(dataSize, 1);
%     clusterCount = 0;
%     
%     for i = 1:dataSize
%         if ~visited(i)
%             visited(i) = true;
%             neighborPts = regionQuery(data, i, epsilon);
%             
%             if numel(neighborPts) < minPts
%                 clusters(i) = -1; % 标记为噪声点
%             else
%                 clusterCount = clusterCount + 1;
%                 expandCluster(data, i, neighborPts, clusterCount, epsilon, minPts, visited, clusters);
%             end
%         end
%     end
% end
%  
% % 密度可达判断
% function neighborPts = regionQuery(data, pointIdx, epsilon)
%     distances = pdist2(data(pointIdx, :), data, 'euclidean');
%     neighborPts = find(distances <= epsilon);
% end
%  
% % 扩展聚类
% function expandCluster(data, pointIdx, neighborPts, clusterCount, epsilon, minPts, visited, clusters)
%     clusters(pointIdx) = clusterCount;
%     i = 1;
%     
%     while i <= numel(neighborPts)
%         currentPt = neighborPts(i);
%         
%         if ~visited(currentPt)
%             visited(currentPt) = true;
%             newNeighborPts = regionQuery(data, currentPt, epsilon);
%             
%             if numel(newNeighborPts) >= minPts
%                 neighborPts = [neighborPts; newNeighborPts'];
%             end
%         end
%         
%         if clusters(currentPt) == 0
%             clusters(currentPt) = clusterCount;
%         end
%         
%         i = i + 1;
%     end
% end

%% 1. 生成随机节点坐标
% 根据粒子的
% 1. 设置随机数生成器的种子
rng(1);

% 2. 生成随机节点坐标
numNodes = 50;  % 总节点数
X = rand(numNodes, 2) * 100;  % 随机生成节点的x, y坐标在 [0,100] 范围内

% 3. 使用 DBSCAN 算法进行聚类 (基于密度的不规则集群)
Eps = 10;  % 邻域半径
MinPts = 3;  % 每个核心点至少的邻居数
[idx, isNoise] = dbscan(X, Eps, MinPts);

% 4. 可视化集群，确保不规则的集群外观
figure;
hold on;
colors = lines(max(idx));  % 使用不同颜色表示不同集群
for i = 1:max(idx)
    clusterNodes = X(idx == i, :);
    if ~isempty(clusterNodes)
        k = boundary(clusterNodes(:,1), clusterNodes(:,2));  % 计算集群的边界
        fill(clusterNodes(k, 1), clusterNodes(k, 2), colors(i,:), 'FaceAlpha', 0.3);  % 绘制不规则边界
        plot(clusterNodes(:, 1), clusterNodes(:, 2), 'o', 'Color', colors(i,:), 'MarkerSize', 8);  % 绘制节点
    end
end

% 5. 绘制噪声点（如果有）
plot(X(isNoise, 1), X(isNoise, 2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);

hold off;
title('基于密度的不规则集群划分');
xlabel('X 坐标');
ylabel('Y 坐标');
grid on;
axis equal;

%% Kmeans
% 1. 生成随机二维数据点
rng(2);
numNodes = 50;  % 数据点的数量
X = rand(numNodes, 2) * 100;  % 随机生成数据点的二维坐标

% 2. 指定聚类的簇数量 k
k = 5;  % 假设希望将数据分为 3 个簇

% 3. 使用 K-means 进行聚类
% idx 是每个点所属的簇的编号, C 是每个簇的质心
[idx, C] = kmeans(X, k);

% 4. 可视化聚类结果
figure;
hold on;
colors = lines(k);  % 定义不同颜色表示不同簇

% 绘制每个簇的点
for i = 1:k
    scatter(X(idx == i, 1), X(idx == i, 2), 100, colors(i,:), 'filled');
end

% 绘制每个簇的质心
plot(C(:, 1), C(:, 2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);

hold off;
title(['K-means 聚类结果, k = ', num2str(k)]);
xlabel('X 坐标');
ylabel('Y 坐标');
grid on;
axis equal;
