% DBSCAN算法
function [clusters, clusterCount] = dbscan(data, epsilon, minPts)
    dataSize = size(data, 1);
    visited = false(dataSize, 1);
    clusters = zeros(dataSize, 1);
    clusterCount = 0;
    
    for i = 1:dataSize
        if ~visited(i)
            visited(i) = true;
            neighborPts = regionQuery(data, i, epsilon);
            
            if numel(neighborPts) < minPts
                clusters(i) = -1; % 标记为噪声点
            else
                clusterCount = clusterCount + 1;
                expandCluster(data, i, neighborPts, clusterCount, epsilon, minPts, visited, clusters);
            end
        end
    end
end
 
% 密度可达判断
function neighborPts = regionQuery(data, pointIdx, epsilon)
    distances = pdist2(data(pointIdx, :), data, 'euclidean');
    neighborPts = find(distances <= epsilon);
end
 
% 扩展聚类
function expandCluster(data, pointIdx, neighborPts, clusterCount, epsilon, minPts, visited, clusters)
    clusters(pointIdx) = clusterCount;
    i = 1;
    
    while i <= numel(neighborPts)
        currentPt = neighborPts(i);
        
        if ~visited(currentPt)
            visited(currentPt) = true;
            newNeighborPts = regionQuery(data, currentPt, epsilon);
            
            if numel(newNeighborPts) >= minPts
                neighborPts = [neighborPts; newNeighborPts'];
            end
        end
        
        if clusters(currentPt) == 0
            clusters(currentPt) = clusterCount;
        end
        
        i = i + 1;
    end
end
 
% 示例数据
data = [
    1, 2;
    1.5, 1.8;
    5, 8;
    8, 8;
    1, 0.6;
    9, 11
];
 
% 设置 DBSCAN 参数
epsilon = 2;
minPts = 2;
 
% 运行 DBSCAN
[clusters, clusterCount] = dbscan(data, epsilon, minPts);
 
% 打印聚类结果
disp('聚类结果：');
disp(clusters);
disp('聚类数目：');
disp(clusterCount);