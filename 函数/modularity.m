%% 模块度
function [rho] = modularity(relationshipMatrix,current_matrix,e)
% 网络边权 柯氏距离 输入邻接矩阵 A (n x n)
A = current_matrix.*e;  % 示例矩阵

% 当节点i和节点j在同一个集群的时候，delta(i,j)等于1，否则等于0 
delta =relationshipMatrix;  % 示例矩阵

% 计算网络中的边数 m
m = sum(A(:)) / 2;

% 计算每个节点的度 k_i
k = sum(A, 2);

% 初始化 rho
rho = 0;

% 计算公式中的 rho
for i = 1:size(A, 1)
    for j = 1:size(A, 1)
        rho = rho + (A(i, j) - (k(i) * k(j)) / (2 * m)) * delta(i,j);
    end
end

rho = (1 / (2 * m)) * rho;

% 输出结果
disp(['The value of ρ is: ', num2str(rho)]);
