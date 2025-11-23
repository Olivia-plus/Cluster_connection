function randomArray = GenerateRandomArray(size, uniqueCount)
    % 输入参数：
    % size: 最终数组的大小（例如 20）
    % uniqueCount: 数组中独特随机数的数量

    % Step 1: 生成 uniqueCount 个随机数
    uniqueNumbers = randi([1, 10], uniqueCount, 1); % 假设随机数范围在 1 到 10 之间
    % Step 2: 重复每个随机数至少 2 次
    repeatedNumbers = repmat(uniqueNumbers, 2, 1); % 重复 2 次

    % Step 3: 随机打乱顺序
    randomArray = repeatedNumbers(randperm(length(repeatedNumbers)));

    % 如果数组大小超过 size，截取前 size 个元素
    if length(randomArray) > size
        randomArray = randomArray(1:size);
    end
    
    % 如果数组大小不足 size，则可以添加更多的随机数
    while length(randomArray) < size
        additionalNumber = randi([1, 10]); % 生成额外的随机数
        randomArray = [randomArray; additionalNumber]; % 添加到数组中
    end

    % 最终随机化顺序
    randomArray = randomArray(randperm(length(randomArray)));
end


% function [particles, clusterIndices] = InitializeParticles(num_buildings, max_num_cluster, x, y)
%     % 生成初始粒子的函数，并返回粒子及其集群划分结果
%     
%     % 生成初始粒子随机数组
%     particles = GenerateRandomArray(num_buildings, max_num_cluster);
%     
%     % 将粒子的坐标与随机数组结合
%     particlesWithCoordinates = [particles', x', y']; % 每行是 [随机数, x, y]
%     
%     % 使用 K-means 算法进行聚类
%     [clusterIndices, ~] = kmeans(particlesWithCoordinates(:, 1:2), max_num_cluster);
%     
%     % 返回粒子和聚类结果
%     particles = particlesWithCoordinates; % 更新粒子为带坐标的粒子
% end
% 
% function randomArray = GenerateRandomArray(size, uniqueCount)
%     % 生成独特随机数的数组
%     uniqueNumbers = randi([1, 10], uniqueCount, 1); % 随机数范围在 1 到 10 之间
%     repeatedNumbers = repmat(uniqueNumbers, 2, 1); % 重复 2 次
%     randomArray = repeatedNumbers(randperm(length(repeatedNumbers))); % 随机打乱顺序
% 
%     % 如果数组大小超过 size，截取前 size 个元素
%     if length(randomArray) > size
%         randomArray = randomArray(1:size);
%     end
%     
%     % 如果数组大小不足 size，则添加更多的随机数
%     while length(randomArray) < size
%         additionalNumber = randi([1, 10]); % 生成额外的随机数
%         randomArray = [randomArray; additionalNumber]; % 添加到数组中
%     end
% 
%     % 最终随机化顺序
%     randomArray = randomArray(randperm(length(randomArray)));
% end







