function process_all_trees(n)
    % 生成所有普吕弗序列，并将其转换为树
    prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
    generate_all_trees(prufer_sequence, 1, n);
end

function generate_all_trees(seq, idx, n)
    % 递归生成普吕弗序列并处理每棵树
    if idx > length(seq)
        tree = prufer_to_tree(seq, n);
        disp('Generated Tree:');
        disp(tree);
        
        % 调用计算权值的函数
        assign_tree_edge_weights_correct(n, tree);
        
        % 计算并显示每棵树的边权和
        edge_sum = calculate_edge_sum(tree);
        disp(['Sum of edge weights: ', num2str(edge_sum)]);
        
    else
        for i = 1:n
            seq(idx) = i;
            generate_all_trees(seq, idx + 1, n);
        end
    end
end

function edge_sum = calculate_edge_sum(tree)
    % 计算每棵树的边权和
    edge_sum = sum(tree(:, 3));  % tree 的第三列表示边的权值
end

function tree = prufer_to_tree(prufer, n)
    % 将普吕弗序列转换为生成树，返回树的边
    degree = ones(1, n);  % 初始化所有节点的度数为1
    for i = 1:length(prufer)
        degree(prufer(i)) = degree(prufer(i)) + 1;  % 序列中的每个节点度数加1
    end

    tree = [];  % 用于存储树的边
    for i = 1:length(prufer)
        for j = 1:n
            if degree(j) == 1
                tree = [tree; j, prufer(i), 0];  % 将节点 j 和序列中的节点连接，初始权值为0
                degree(j) = degree(j) - 1;
                degree(prufer(i)) = degree(prufer(i)) - 1;
                break;
            end
        end
    end

    % 最后剩下两个节点，直接相连
    remaining = find(degree == 1);
    tree = [tree; remaining(1), remaining(2), 0];  % 初始权值为0
end

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
    edge_weights(end) = max(node_values(remaining_edge(1:2))); % 取两侧中较大的值
    
    % 将权值添加到 tree 的第三列
    tree(:, 3) = edge_weights;
    
    % 输出最终的边和对应的权值
    disp('Edges with Weights:');
    disp(tree);
end

function [updated_tree, updated_node_values, updated_edge_weights] = process_leaves_correct(tree, node_values, edge_weights)
    % 找到度为 1 的节点（叶子节点）
    degree = zeros(max(max(tree(:, 1:2))), 1);  % 初始化所有节点的度数
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


