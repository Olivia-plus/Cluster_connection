function [e] = complementarity(net_load,num_buildings)
% 计算互补矩阵
% 假设有 m 个节点，T 个时间点
% m = 5; % 节点数量
% T = 10; % 时间点数量

% 随机生成 P 矩阵 (m x T), 代表每个节点在每个时间点的净功率
% 你可以替换为实际的净功率数据
% P = rand(m, T); 

% 初始化互补度边权矩阵
e = zeros(num_buildings, num_buildings);

% 计算每个节点对 (i, j) 的边权
for i = 1:num_buildings
    for j = 1:num_buildings
        if i ~= j
            % 计算分子
            numerator = sum(net_load{i})+ sum(net_load{j});
            % 计算所有节点对中的最大值（用于分母）
            max_value = 0;
            for ii = 1:num_buildings
                for jj = ii+1:num_buildings
                    current_value = sum(net_load{ii}) + sum(net_load{jj});
                    if abs(current_value) > max_value
                        max_value = abs(current_value);
                    end
                end
            end
            % 计算 e_ij
            e(i, j) = 1 - abs(numerator) / max_value;
        end
    end
end

% 输出边权矩阵
disp(e);

end

