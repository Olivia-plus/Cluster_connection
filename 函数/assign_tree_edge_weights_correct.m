function assign_tree_edge_weights_correct(n, tree)
    % 初始化每个节点的值为它的编号
    node_values = 1:n;
    
    % 初始化每条边的权值
    edge_weights = zeros(size(tree, 1), 1);
    
    % 逐层处理叶子节点，递归赋值边的权值
    while size(tree, 1) > 1
        [tree, node_values, edge_weights] = process_leaves_correct(tree, node_values, edge_weights);
    end
    
    % 处理最后一条边
    remaining_edge = tree(1, :);
    edge_weights(end) = max(node_values(remaining_edge)); % 取两侧中较大的值
    
    % 输出最终的边和对应的权值
    disp('Edges:');
    disp(tree);
    disp('Edge Weights:');
    disp(edge_weights);
end

function [updated_tree, updated_node_values, updated_edge_weights] = process_leaves_correct(tree, node_values, edge_weights)
    % 找到度为 1 的节点（叶子节点）
    degree = zeros(max(max(tree)), 1);  % 初始化所有节点的度数
    for i = 1:size(tree, 1)
        degree(tree(i, 1)) = degree(tree(i, 1)) + 1;
        degree(tree(i, 2)) = degree(tree(i, 2)) + 1;
    end
    leaf_nodes = find(degree == 1);  % 度为 1 的节点即叶子节点
    
    % 处理每个叶子节点
    for i = 1:length(leaf_nodes)
        leaf = leaf_nodes(i);
        
        % 找到这个叶子节点连接的边
        for j = 1:size(tree, 1)
            if tree(j, 1) == leaf || tree(j, 2) == leaf
                % 获取连接的另一个节点
                connected_node = tree(j, 1) + tree(j, 2) - leaf;
                
                % 如果这条边没有赋过权值，给它赋值
                if edge_weights(j) == 0
                    edge_weights(j) = node_values(leaf);  % 当前叶子的值赋给该边
                end
                
                % 删除这个叶子节点和边
                tree(j, :) = [];
                break;
            end
        end
    end
    
    updated_tree = tree;
    updated_node_values = node_values;
    updated_edge_weights = edge_weights;
end

% 示例：给树 1-2-3-4 赋值
n = 4;
tree = [1 2; 2 3; 3 4];  % 树结构，表示边
assign_tree_edge_weights_correct(n, tree);
