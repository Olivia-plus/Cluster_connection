function find_all_tree_weights(n)
    % 生成所有普吕弗序列，并针对每棵树进行遍历
    all_prufer_sequences = generate_all_prufer_sequences(n);
    
    % 遍历每个普吕弗序列生成的树
    for i = 1:size(all_prufer_sequences, 1)
        prufer_seq = all_prufer_sequences(i, :);
        tree = prufer_to_tree(prufer_seq, n);
        disp('Tree structure:');
        disp(tree);
        
        % 初始化所有边的权值
        edge_weights = zeros(size(tree, 1), 1);
        
        % 在这棵树中按顺序确定每条边的权值
        current_tree = tree;
        remaining_nodes = n;
        while remaining_nodes > 1
            % 找到当前树中的度为1的节点
            [current_tree, edge_weights] = assign_weights_to_leaves(current_tree, edge_weights);
            remaining_nodes = size(current_tree, 1) + 1;  % 更新剩余节点数
        end
        
        disp('Edge weights:');
        disp(edge_weights);
    end
end

function all_sequences = generate_all_prufer_sequences(n)
    % 生成所有长度为 n-2 的普吕弗序列
    sequence_length = n - 2;
    all_sequences = perms(1:n);  % 所有节点的排列
    all_sequences = all_sequences(:, 1:sequence_length);  % 取前 n-2 项
end

function tree = prufer_to_tree(prufer, n)
    % 将普吕弗序列转换为树的边列表
    degree = ones(1, n);  % 所有节点的初始度数为1
    for i = 1:length(prufer)
        degree(prufer(i)) = degree(prufer(i)) + 1;
    end
    
    tree = [];  % 初始化边列表
    for i = 1:length(prufer)
        for j = 1:n
            if degree(j) == 1
                tree = [tree; j, prufer(i)];  % 连接叶子节点与序列节点
                degree(j) = degree(j) - 1;
                degree(prufer(i)) = degree(prufer(i)) - 1;
                break;
            end
        end
    end
    
    % 剩下最后两个节点，直接连接
    remaining = find(degree == 1);
    tree = [tree; remaining(1), remaining(2)];
end

function [updated_tree, edge_weights] = assign_weights_to_leaves(tree, edge_weights)
    % 找到树中度为1的节点并分配权值
    degree = zeros(max(max(tree)), 1);  % 初始化节点度数
    for i = 1:size(tree, 1)
        degree(tree(i, 1)) = degree(tree(i, 1)) + 1;
        degree(tree(i, 2)) = degree(tree(i, 2)) + 1;
    end
    
    % 找到度为1的节点
    leaf_nodes = find(degree == 1);
    
    % 分配权值给这些叶子节点连接的边
    for i = 1:size(tree, 1)
        if ismember(tree(i, 1), leaf_nodes) || ismember(tree(i, 2), leaf_nodes)
            % 假设当前边的权值为叶子节点的序号
            edge_weights(i) = min(tree(i, :));
        end
    end
    
    % 删除这些叶子节点和相连的边
    updated_tree = tree;
    for i = 1:length(leaf_nodes)
        updated_tree(any(updated_tree == leaf_nodes(i), 2), :) = [];
    end
end


