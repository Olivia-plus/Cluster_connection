% %% 参数设置
% numPoints = 20;          % 随机生成的点数
% numClusters = 5;         % 期望的集群数量
% numParticles = 30;       % 粒子群中的粒子数量
% maxIter = 100;           % 最大迭代次数
% w = 0.5;                 % 惯性权重
% c1 = 1.5;                % 个体加速系数
% c2 = 1.5;                % 群体加速系数
% velocityLimit = 0.1;     % 粒子的速度限制
% 
% %% 随机生成数据点
% X = rand(numPoints, 2) * 10;  % 随机生成 numPoints 个二维点，范围 [0, 10]
% 
% % 显示初始数据点分布
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% title('随机生成的节点分布');
% xlabel('X');
% ylabel('Y');
% hold on;
% 
% %% 初始化粒子群
% particles = cell(numParticles, 1);  % 每个粒子代表一个质心集
% velocities = cell(numParticles, 1); % 每个粒子的速度
% 
% % 初始化每个粒子的质心和速度
% for i = 1:numParticles
%     particles{i} = X(randperm(numPoints, numClusters), :);  % 随机选择初始质心
%     velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
% end
% 
% % 初始化最优解
% pbest = particles;       % 个体最优解
% pbestFitness = inf(numParticles, 1); % 个体最优适应度
% gbest = particles{1};    % 全局最优解
% gbestFitness = inf;      % 全局最优适应度
% maxIter=100;
% %% 迭代更新粒子群
% for iter = 1:maxIter
%     for p = 1:numParticles
%         % 获取当前粒子的质心
%         C = particles{p};
% 
%         % 计算每个点到质心的距离
%         distances = pdist2(X, C);
%         [~, clusterIdx] = min(distances, [], 2);  % 分配每个点到最近的质心
%         
%         % 检查当前的集群划分是否满足条件
%         valid = checkClusterValidity(X, clusterIdx, numClusters);
%         
%         % 计算适应度 (集群内的距离总和作为适应度)
%         if valid
%             fitness = computeFitness(X, C, clusterIdx);
%         else
%             fitness = inf;  % 不合法的解给高适应度值，避免选择
%         end
%         
%         % 更新个体最优
%         if fitness < pbestFitness(p)
%             pbest{p} = particles{p};
%             pbestFitness(p) = fitness;
%         end
%         
%         % 更新全局最优
%         if fitness < gbestFitness
%             gbest = particles{p};
%             gbestFitness = fitness;
%         end
%     end
%     
%     % 更新粒子的位置和速度
%     for p = 1:numParticles
%         velocities{p} = w * velocities{p} ...
%             + c1 * rand * (pbest{p} - particles{p}) ...
%             + c2 * rand * (gbest - particles{p});
%         particles{p} = particles{p} + velocities{p};  % 更新粒子的位置
%     end
%     
%     % 每 10 次迭代可视化一次当前最优解
%     if mod(iter, 10) == 0
%         figure;
%         scatter(X(:,1), X(:,2), 'filled');
%         hold on;
%         scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
%         title(['迭代次数: ', num2str(iter), ' 最优适应度: ', num2str(gbestFitness)]);
%         xlabel('X');
%         ylabel('Y');
%         hold off;
%         pause(0.5);
%     end
%     
%     % 检查是否达到收敛条件
%     if gbestFitness < 1e-5
%         break;
%     end
% end
% 
% % 输出最终最优解
% disp('最终最优集群划分质心:');
% disp(gbest);
% 
% % 显示最终结果
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% hold on;
% scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
% title('最终最优集群划分');
% xlabel('X');
% ylabel('Y');
% function valid = checkClusterValidity(X, clusterIdx, numClusters)
%     % 检查每个集群是否至少有两个点
%     for i = 1:numClusters
%         if sum(clusterIdx == i) < 2
%             valid = false;
%             return;
%         end
%     end
%     
%     % 检查集群之间是否重叠，避免质心过近
%     clusterCenters = zeros(numClusters, 2);
%     for i = 1:numClusters
%         clusterCenters(i, :) = mean(X(clusterIdx == i, :), 1);
%     end
%     
%     % 计算质心之间的距离
%     if min(pdist(clusterCenters)) < 0.5  % 集群之间的最小距离，阈值可调整
%         valid = false;
%         return;
%     end
%     
%     valid = true;
% end
% function fitness = computeFitness(X, C, clusterIdx)
%     fitness = 0;
%     numClusters = size(C, 1);
%     for i = 1:numClusters
%         % 计算每个集群内点到质心的距离和
%         pointsInCluster = X(clusterIdx == i, :);
%         if ~isempty(pointsInCluster)
%             fitness = fitness + sum(sum((pointsInCluster - C(i, :)).^2));
%         end
%     end
% end




% %% 参数设置
% rng(3)
% numPoints = 50;          % 随机生成的点数
% numClusters = 5;         % 集群数量
% numParticles = 30;       % 粒子群数量
% maxIter = 100;           % 最大迭代次数
% w = 0.8;                 % 惯性权重
% c1 = 1.5;                % 个体加速系数
% c2 = 1.5;                % 群体加速系数
% velocityLimit = 4;     % 粒子速度限制
% 
% %% 随机生成数据点
% X = rand(numPoints, 2) * 10;  % 随机生成 numPoints 个二维点，范围 [0, 10]
% 
% % 显示初始点分布
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% title('随机生成的节点分布');
% xlabel('X');
% ylabel('Y');
% hold on;
% 
% %% 初始化粒子群
% particles = cell(numParticles, 1);  % 每个粒子表示不同的质心集
% velocities = cell(numParticles, 1); % 每个粒子的速度
% 
% for i = 1:numParticles
%     particles{i} = X(randperm(numPoints, numClusters), :);  % 随机选择质心
%     velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
% end
% 
% % 个体最优解
% pbest = particles;       % 个体最优解
% pbestFitness = inf(numParticles, 1); % 个体最优适应度
% gbest = particles{1};    % 全局最优解
% gbestFitness = inf;      % 全局最优适应度
% 
% %% 迭代更新质心位置
% for iter = 1:maxIter
%     for p = 1:numParticles
%         % 获取当前粒子的质心
%         C = particles{p};
%         
%         % 计算每个点到质心的距离
%         distances = pdist2(X, C);
%         [~, clusterIdx] = min(distances, [], 2);  % 将每个点分配到最近的质心
%         
%         % 检查集群是否满足要求：没有孤立节点，集群不重叠
%         valid = checkClusterValidity(X, clusterIdx, numClusters);
%         
%         % 计算适应度（集群内距离之和）
%         if valid
%             fitness = computeFitness(X, C, clusterIdx);
%         else
%             fitness = inf;  % 不合法的集群划分给予惩罚
%         end
%         
%         % 更新个体最优解
%         if fitness < pbestFitness(p)
%             pbest{p} = particles{p};
%             pbestFitness(p) = fitness;
%         end
%         
%         % 更新全局最优解
%         if fitness < gbestFitness
%             gbest = particles{p};
%             gbestFitness = fitness;
%         end
%     end
%     
%     % 更新质心位置
%     for p = 1:numParticles
%         velocities{p} = w * velocities{p} ...
%             + c1 * rand * (pbest{p} - particles{p}) ...
%             + c2 * rand * (gbest - particles{p});  % 粒子速度更新
%         particles{p} = particles{p} + velocities{p};  % 更新粒子质心位置
%     end
%     
%     % 可视化当前最优集群划分
%     if mod(iter, 10) == 0
%         figure;
%         scatter(X(:,1), X(:,2), 'filled');
%         hold on;
%         scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
%         title(['迭代次数: ', num2str(iter), ' 最优适应度: ', num2str(gbestFitness)]);
%         xlabel('X');
%         ylabel('Y');
%         hold off;
%         pause(0.5);
%     end
%     
%     % 检查是否达到收敛条件
%     if gbestFitness < 1e-5
%         break;
%     end
% end
% 
% % 显示最终集群划分结果
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% hold on;
% scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
% title('最终最优集群划分');
% xlabel('X');
% ylabel('Y');
% 
% %% 辅助函数：检查集群有效性
% function valid = checkClusterValidity(X, clusterIdx, numClusters)
%     for i = 1:numClusters
%         if sum(clusterIdx == i) < 2
%             valid = false;
%             return;
%         end
%     end
%     
%     clusterCenters = zeros(numClusters, 2);
%     for i = 1:numClusters
%         clusterCenters(i, :) = mean(X(clusterIdx == i, :), 1);
%     end
%     
%     if min(pdist(clusterCenters)) < 0.5
%         valid = false;
%         return;
%     end
%     valid = true;
% end
% 
% %% 辅助函数：计算适应度
% function fitness = computeFitness(X, C, clusterIdx)
%     fitness = 0;
%     numClusters = size(C, 1);
%     for i = 1:numClusters
%         pointsInCluster = X(clusterIdx == i, :);
%         if ~isempty(pointsInCluster)
%             fitness = fitness + sum(sum((pointsInCluster - C(i, :)).^2));
%         end
%     end
% end
% 
% %% 参数设置
% rng(4)
% numPoints = 50;          % 随机生成的点数
% numClusters = 5;         % 集群数量
% numParticles = 30;       % 粒子群数量
% maxIter = 30;           % 最大迭代次数
% w = 0.5;                 % 惯性权重
% c1 = 1.5;                % 个体加速系数
% c2 = 1.5;                % 群体加速系数
% velocityLimit = 1;     % 粒子速度限制
% 
% %% 随机生成数据点
% X = rand(numPoints, 2) * 10;  % 随机生成 numPoints 个二维点，范围 [0, 10]
% 
% % 显示初始点分布
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% title('随机生成的节点分布');
% xlabel('X');
% ylabel('Y');
% hold on;
% 
% %% 初始化粒子群
% particles = cell(numParticles, 1);  % 每个粒子表示不同的质心集
% velocities = cell(numParticles, 1); % 每个粒子的速度
% 
% for i = 1:numParticles
%     particles{i} = X(randperm(numPoints, numClusters), :);  % 随机选择质心
%     velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
% end
% 
% % 个体最优解
% pbest = particles;       % 个体最优解
% pbestFitness = inf(numParticles, 1); % 个体最优适应度
% gbest = particles{1};    % 全局最优解
% gbestFitness = inf;      % 全局最优适应度
% 
% % 记录适应度变化
% fitnessHistory = zeros(maxIter, 1);  % 适应度变化记录
% 
% %% 迭代更新质心位置
% for iter = 1:maxIter
%     for p = 1:numParticles
%         % 获取当前粒子的质心
%         C = particles{p};
%         
%         % 计算每个点到质心的距离
%         distances = pdist2(X, C);
%         [~, clusterIdx] = min(distances, [], 2);  % 将每个点分配到最近的质心
%         
%         % 检查集群是否满足要求：没有孤立节点，集群不重叠
%         valid = checkClusterValidity(X, clusterIdx, numClusters);
%         
%         % 计算适应度（集群内距离之和）
%         if valid
%             fitness = computeFitness(X, C, clusterIdx);
%         else
%             fitness = inf;  % 不合法的集群划分给予惩罚
%         end
%         
%         % 更新个体最优解
%         if fitness < pbestFitness(p)
%             pbest{p} = particles{p};
%             pbestFitness(p) = fitness;
%         end
%         
%         % 更新全局最优解
%         if fitness < gbestFitness
%             gbest = particles{p};
%             gbestFitness = fitness;
%         end
%     end
%     
%     % 记录当前迭代的适应度
%     fitnessHistory(iter) = gbestFitness;
%     
%     % 更新质心位置
%     for p = 1:numParticles
%         velocities{p} = w * velocities{p} ...
%             + c1 * rand * (pbest{p} - particles{p}) ...
%             + c2 * rand * (gbest - particles{p});  % 粒子速度更新
%         particles{p} = particles{p} + velocities{p};  % 更新粒子质心位置
%     end
%     
%     % 可视化当前最优集群划分
%     if mod(iter, 10) == 0
%         figure;
%         scatter(X(:,1), X(:,2), 'filled');
%         hold on;
%         scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
%         title(['迭代次数: ', num2str(iter), ' 最优适应度: ', num2str(gbestFitness)]);
%         xlabel('X');
%         ylabel('Y');
%         hold off;
%         pause(0.5);
%     end
%     
%     % 检查是否达到收敛条件
%     if gbestFitness < 1e-5
%         break;
%     end
% end
% 
% % 显示最终集群划分结果
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% hold on;
% scatter(gbest(:,1), gbest(:,2), 100, 'rx', 'LineWidth', 2);
% title('最终最优集群划分');
% xlabel('X');
% ylabel('Y');
% 
% % 绘制适应度变化图
% figure;
% plot(fitnessHistory(1:iter), 'LineWidth', 2);
% title('适应度变化图');
% xlabel('迭代次数');
% ylabel('适应度');
% grid on;
% 
% %% 辅助函数：检查集群有效性
% function valid = checkClusterValidity(X, clusterIdx, numClusters)
%     for i = 1:numClusters
%         if sum(clusterIdx == i) < 2
%             valid = false;
%             return;
%         end
%     end
%     
%     clusterCenters = zeros(numClusters, 2);
%     for i = 1:numClusters
%         clusterCenters(i, :) = mean(X(clusterIdx == i, :), 1);
%     end
%     
%     if min(pdist(clusterCenters)) < 0.5
%         valid = false;
%         return;
%     end
%     valid = true;
% end
% 
% %% 辅助函数：计算适应度
% function fitness = computeFitness(X, C, clusterIdx)
%     fitness = 0;
%     numClusters = size(C, 1);
%     for i = 1:numClusters
%         pointsInCluster = X(clusterIdx == i, :);
%         if ~isempty(pointsInCluster)
%             fitness = fitness + sum(sum((pointsInCluster - C(i, :)).^2));
%         end
%     end
% end
% %% 参数设置
% numPoints = 50;          % 随机生成的点数
% numClusters = 5;         % 集群数量
% numParticles = 30;       % 粒子群数量
% maxIter = 100;           % 最大迭代次数
% w = 0.5;                 % 惯性权重
% c1 = 1.5;                % 个体加速系数
% c2 = 1.5;                % 群体加速系数
% velocityLimit = 0.1;     % 粒子速度限制
% 
% %% 随机生成数据点
% X = rand(numPoints, 2) * 10;  % 随机生成 numPoints 个二维点，范围 [0, 10]
% 
% % 显示初始点分布
% figure;
% scatter(X(:,1), X(:,2), 'filled');
% title('随机生成的节点分布');
% xlabel('X');
% ylabel('Y');
% hold on;
% 
% %% 初始化粒子群
% particles = cell(numParticles, 1);  % 每个粒子表示不同的质心集
% velocities = cell(numParticles, 1); % 每个粒子的速度
% 
% for i = 1:numParticles
%     particles{i} = X(randperm(numPoints, numClusters), :);  % 随机选择质心
%     velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
% end
% 
% % 个体最优解
% pbest = particles;       % 个体最优解
% pbestFitness = inf(numParticles, 1); % 个体最优适应度
% gbest = particles{1};    % 全局最优解
% gbestFitness = inf;      % 全局最优适应度
% 
% % 记录适应度变化
% fitnessHistory = zeros(maxIter, 1);  % 适应度变化记录
% 
% %% 迭代更新质心位置
% for iter = 1:maxIter
%     for p = 1:numParticles
%         % 获取当前粒子的质心
%         C = particles{p};
%         
%         % 计算每个点到质心的距离
%         distances = pdist2(X, C);
%         [~, clusterIdx] = min(distances, [], 2);  % 将每个点分配到最近的质心
%         
%         % 检查集群是否满足要求：没有孤立节点，集群不重叠
%         valid = checkClusterValidity(X, clusterIdx, numClusters);
%         
%         % 计算适应度（集群内距离之和）
%         if valid
%             fitness = computeFitness(X, C, clusterIdx);
%         else
%             fitness = inf;  % 不合法的集群划分给予惩罚
%         end
%         
%         % 更新个体最优解
%         if fitness < pbestFitness(p)
%             pbest{p} = particles{p};
%             pbestFitness(p) = fitness;
%         end
%         
%         % 更新全局最优解
%         if fitness < gbestFitness
%             gbest = particles{p};
%             gbestFitness = fitness;
%         end
%     end
%     
%     % 记录当前迭代的适应度
%     fitnessHistory(iter) = gbestFitness;
%     
%     % 更新质心位置
%     for p = 1:numParticles
%         velocities{p} = w * velocities{p} ...
%             + c1 * rand * (pbest{p} - particles{p}) ...
%             + c2 * rand * (gbest - particles{p});  % 粒子速度更新
%         particles{p} = particles{p} + velocities{p};  % 更新粒子质心位置
%     end
%     
%     % 可视化当前最优集群划分（每隔10次迭代）
%     if mod(iter, 10) == 0
%         visualizeClusters(X, gbest, clusterIdx, numClusters);
%         title(['迭代次数: ', num2str(iter), ' 最优适应度: ', num2str(gbestFitness)]);
%         pause(0.5);
%     end
%     
%     % 检查是否达到收敛条件
%     if gbestFitness < 1e-5
%         break;
%     end
% end
% 
% % 显示最终集群划分结果
% visualizeClusters(X, gbest, clusterIdx, numClusters);
% title('最终最优集群划分');
% 
% % 绘制适应度变化图
% figure;
% plot(fitnessHistory(1:iter), 'LineWidth', 2);
% title('适应度变化图');
% xlabel('迭代次数');
% ylabel('适应度');
% grid on;
% 
% %% 辅助函数：检查集群有效性
% function valid = checkClusterValidity(X, clusterIdx, numClusters)
%     for i = 1:numClusters
%         if sum(clusterIdx == i) < 2
%             valid = false;
%             return;
%         end
%     end
%     
%     clusterCenters = zeros(numClusters, 2);
%     for i = 1:numClusters
%         clusterCenters(i, :) = mean(X(clusterIdx == i, :), 1);
%     end
%     
%     if min(pdist(clusterCenters)) < 0.5
%         valid = false;
%         return;
%     end
%     valid = true;
% end
% 
% %% 辅助函数：计算适应度
% function fitness = computeFitness(X, C, clusterIdx)
%     fitness = 0;
%     numClusters = size(C, 1);
%     for i = 1:numClusters
%         pointsInCluster = X(clusterIdx == i, :);
%         if ~isempty(pointsInCluster)
%             fitness = fitness + sum(sum((pointsInCluster - C(i, :)).^2));
%         end
%     end
% end
% 
% %% 辅助函数：可视化集群划分
% function visualizeClusters(X, C, clusterIdx, numClusters)
%     figure;
%     colors = lines(numClusters);  % 使用不同的颜色标识不同集群
%     hold on;
%     
%     for i = 1:numClusters
%         clusterPoints = X(clusterIdx == i, :);  % 选出属于第 i 个集群的点
%         scatter(clusterPoints(:, 1), clusterPoints(:, 2), 50, colors(i, :), 'filled');
%         
%         % 画出集群的边界（凸包）
%         if size(clusterPoints, 1) > 2
%             k = convhull(clusterPoints(:, 1), clusterPoints(:, 2));
%             plot(clusterPoints(k, 1), clusterPoints(k, 2), 'Color', colors(i, :), 'LineWidth', 2);
%         end
%         
%         scatter(C(i, 1), C(i, 2), 100, 'x', 'LineWidth', 2, 'MarkerEdgeColor', colors(i, :));  % 质心
%     end
%     
%     title('集群划分结果');
%     xlabel('X');
%     ylabel('Y');
%     hold off;
% end
%% 参数设置
numPoints = 50;          % 随机生成的点数
numClusters = 8;         % 集群数量
numParticles = 10;       % 粒子群数量
maxIter = 25;           % 最大迭代次数
w = 0.5;                 % 惯性权重
c1 = 1.5;                % 个体加速系数
c2 = 1.5;                % 群体加速系数
velocityLimit = 0.1;     % 粒子速度限制

%% 随机生成数据点
X = rand(numPoints, 2) * 10;  % 随机生成 numPoints 个二维点，范围 [0, 10]

% 显示初始点分布
figure;
scatter(X(:,1), X(:,2), 'filled');
title('随机生成的节点分布');
xlabel('X');
ylabel('Y');
hold on;

%% 初始化粒子群
particles = cell(numParticles, 1);  % 每个粒子表示不同的质心集
velocities = cell(numParticles, 1); % 每个粒子的速度

for i = 1:numParticles
    particles{i} = X(randperm(numPoints, numClusters), :);  % 随机选择质心
    velocities{i} = randn(numClusters, 2) * velocityLimit;  % 初始化速度
end

% 个体最优解
pbest = particles;       % 个体最优解
pbestFitness = inf(numParticles, 1); % 个体最优适应度
gbest = particles{1};    % 全局最优解
gbestFitness = inf;      % 全局最优适应度

% 记录适应度变化
fitnessHistory = zeros(maxIter, 1);  % 适应度变化记录

%% 迭代更新质心位置
for iter = 1:maxIter
    for p = 1:numParticles
        % 获取当前粒子的质心
        C = particles{p};
        
        % 计算每个点到质心的距离
        distances = pdist2(X, C);
        [~, clusterIdx] = min(distances, [], 2);  % 将每个点分配到最近的质心
        
        % 检查集群是否满足要求：没有孤立节点，集群不重叠
        valid = checkClusterValidity(X, clusterIdx, numClusters);
        
        % 计算适应度（集群内距离之和）
        if valid
            fitness = computeFitness(X, C, clusterIdx);
        else
            fitness = inf;  % 不合法的集群划分给予惩罚
        end
        
        % 更新个体最优解
        if fitness < pbestFitness(p)
            pbest{p} = particles{p};
            pbestFitness(p) = fitness;
        end
        
        % 更新全局最优解
        if fitness < gbestFitness
            gbest = particles{p};
            gbestFitness = fitness;
        end
    end
    
    % 记录当前迭代的适应度
    fitnessHistory(iter) = gbestFitness;
    
    % 更新质心位置
    for p = 1:numParticles
        velocities{p} = w * velocities{p} ...
            + c1 * rand * (pbest{p} - particles{p}) ...
            + c2 * rand * (gbest - particles{p});  % 粒子速度更新
        particles{p} = particles{p} + velocities{p};  % 更新粒子质心位置
    end
    
    % 可视化第一个粒子的集群划分
    figure(1);
    subplot(1, 2, 1);
    visualizeClusters(X, particles{1}, clusterIdx, numClusters);
    title(['第一个粒子的集群划分 (迭代: ', num2str(iter), ')']);
    
    % 可视化当前最优集群划分（每隔10次迭代）
    subplot(1, 2, 2);
    visualizeClusters(X, gbest, clusterIdx, numClusters);
    title(['全局最优集群划分 (迭代: ', num2str(iter), ')']);
    
    pause(0.5);
    
    % 检查是否达到收敛条件
    if gbestFitness < 1e-5
        break;
    end
end

% 显示最终集群划分结果
figure(2);
visualizeClusters(X, gbest, clusterIdx, numClusters);
title('最终最优集群划分');

% 绘制适应度变化图
figure(3);
plot(fitnessHistory(1:iter), 'LineWidth', 2);
title('适应度变化图');
xlabel('迭代次数');
ylabel('适应度');
grid on;

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
    
    if min(pdist(clusterCenters)) < 0.5
        valid = false;
        return;
    end
    valid = true;
end

%% 辅助函数：计算适应度
function fitness = computeFitness(X, C, clusterIdx)
    fitness = 0;
    numClusters = size(C, 1);
    for i = 1:numClusters
        pointsInCluster = X(clusterIdx == i, :);
        if ~isempty(pointsInCluster)
            fitness = fitness + sum(sum((pointsInCluster - C(i, :)).^2));
        end
    end
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
