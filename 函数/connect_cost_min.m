% 构造函数：传入集群的建筑数据，需要传入三维的建筑功率传输矩阵，以及建筑的数量，建筑的x,y的值
function [adjacency_matrix,min_cost]=connect_cost_min(P,n,x_cluster,y_cluster)
% n = 4;
process_all_trees(n,x_cluster,y_cluster);

% % 示例：生成 n = 4 的所有树并计算每棵树的边权和
% % function process_all_trees(n)
% %     % 生成所有普吕弗序列，并将其转换为树
% %     prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
% %     generate_all_trees(prufer_sequence, 1, n);
% % end
% % 
% % function generate_all_trees(seq, idx, n)
% %     % 递归生成普吕弗序列并处理每棵树
% %     if idx > length(seq)
% %         tree = prufer_to_tree(seq, n);
% %         disp('Generated Tree:');
% %         disp(tree);
% %         
% %         % 调用计算权值的函数
% %         edge_weights = assign_tree_edge_weights_correct(n, tree);
% %         
% %         % 计算并显示每棵树的边权和
% %         edge_sum = sum(edge_weights);  % 计算边权和
% %         disp(['Sum of edge weights: ', num2str(edge_sum)]);
% %         
% %     else
% %         for i = 1:n
% %             seq(idx) = i;
% %             generate_all_trees(seq, idx + 1, n);
% %         end
% %     end
% % end
% % 
% % function tree = prufer_to_tree(prufer, n)
% %     % 将普吕弗序列转换为生成树，返回树的边
% %     degree = ones(1, n);  % 初始化所有节点的度数为1
% %     for i = 1:length(prufer)
% %         degree(prufer(i)) = degree(prufer(i)) + 1;  % 序列中的每个节点度数加1
% %     end
% % 
% %     tree = [];  % 用于存储树的边
% %     for i = 1:length(prufer)
% %         for j = 1:n
% %             if degree(j) == 1
% %                 tree = [tree; j, prufer(i)];  % 将节点 j 和序列中的节点连接
% %                 degree(j) = degree(j) - 1;
% %                 degree(prufer(i)) = degree(prufer(i)) - 1;
% %                 break;
% %             end
% %         end
% %     end
% % 
% %     % 最后剩下两个节点，直接相连
% %     remaining = find(degree == 1);
% %     tree = [tree; remaining(1), remaining(2)];
% % end
% % 
% % function edge_weights = assign_tree_edge_weights_correct(n, tree)
% %     % 初始化每个节点的值为它的编号
% %     node_values = 1:n;
% %     
% %     % 初始化每条边的权值
% %     edge_weights = zeros(size(tree, 1), 1);
% %     
% %     % 逐层处理叶子节点，递归赋值边的权值
% %     while size(tree, 1) > 1
% %         [tree, node_values, edge_weights] = process_leaves_correct(tree, node_values, edge_weights);
% %     end
% %     
% %     % 如果还剩一条边，处理最后一条边
% %     if size(tree, 1) == 1
% %         remaining_edge = tree(1, :);
% %         edge_weights(end) = max(node_values(remaining_edge)); % 取两侧中较大的值
% %     end
% % end
% % 
% % function [updated_tree, updated_node_values, updated_edge_weights] = process_leaves_correct(tree, node_values, edge_weights)
% %     % 找到度为 1 的节点（叶子节点）
% %     degree = zeros(max(max(tree(:, 1:2))), 1);  % 初始化所有节点的度数
% %     for i = 1:size(tree, 1)
% %         degree(tree(i, 1)) = degree(tree(i, 1)) + 1;
% %         degree(tree(i, 2)) = degree(tree(i, 2)) + 1;
% %     end
% %     leaf_nodes = find(degree == 1);  % 度为 1 的节点即叶子节点
% %     
% %     % 处理每个叶子节点
% %     for i = 1:length(leaf_nodes)
% %         leaf = leaf_nodes(i);
% %         
% %         % 找到这个叶子节点连接的边
% %         for j = 1:size(tree, 1)
% %             if tree(j, 1) == leaf || tree(j, 2) == leaf
% %                 % 获取连接的另一个节点
% %                 connected_node = tree(j, 1) + tree(j, 2) - leaf;
% %                 
% %                 % 如果这条边没有赋过权值，给它赋值
% %                 if edge_weights(j) == 0
% %                     edge_weights(j) = node_values(leaf);  % 当前叶子的值赋给该边
% %                 end
% %                 
% %                 % 更新节点值：将叶子节点值与相连节点的值累加，更新相连节点
% %                 node_values(connected_node) = node_values(connected_node) + node_values(leaf);
% %                 
% %                 % 删除这个叶子节点和边
% %                 tree(j, :) = [];
% %                 break;
% %             end
% %         end
% %     end
% %     
% %     updated_tree = tree;
% %     updated_node_values = node_values;
% %     updated_edge_weights = edge_weights;
% % end
% 
% % % 改了输出格式
% % function process_all_trees(n)
% %     % 生成所有普吕弗序列，并将其转换为树
% %     prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
% %     generate_all_trees(prufer_sequence, 1, n);
% % end
% % 
% % function generate_all_trees(seq, idx, n)
% %     % 递归生成普吕弗序列并处理每棵树
% %     if idx > length(seq)
% %         tree = prufer_to_tree(seq, n);
% %         disp('Generated Tree:');
% %         disp(tree);  % 显示生成的树（边列表）
% % 
% %         % 调用计算权值的函数
% %         edge_weights = assign_tree_edge_weights_correct(n, tree);
% %         
% %         % 显示树的边及权重
% %         disp('Edges with Weights:');
% %         for i = 1:size(tree, 1)
% %             fprintf('%5d %5d %5d\n', tree(i, 1), tree(i, 2), edge_weights(i));
% %         end
% %         
% %         % 计算并显示每棵树的边权和
% %         edge_sum = sum(edge_weights);  % 计算边权和
% %         disp(['Sum of edge weights: ', num2str(edge_sum)]);
% %         disp('----------------------------');
% %         
% %     else
% %         for i = 1:n
% %             seq(idx) = i;
% %             generate_all_trees(seq, idx + 1, n);
% %         end
% %     end
% % end
% % 
% % function tree = prufer_to_tree(prufer, n)
% %     % 将普吕弗序列转换为生成树，返回树的边
% %     degree = ones(1, n);  % 初始化所有节点的度数为1
% %     for i = 1:length(prufer)
% %         degree(prufer(i)) = degree(prufer(i)) + 1;  % 序列中的每个节点度数加1
% %     end
% % 
% %     tree = [];  % 用于存储树的边
% %     for i = 1:length(prufer)
% %         for j = 1:n
% %             if degree(j) == 1
% %                 tree = [tree; j, prufer(i)];  % 将节点 j 和序列中的节点连接
% %                 degree(j) = degree(j) - 1;
% %                 degree(prufer(i)) = degree(prufer(i)) - 1;
% %                 break;
% %             end
% %         end
% %     end
% % 
% %     % 最后剩下两个节点，直接相连
% %     remaining = find(degree == 1);
% %     tree = [tree; remaining(1), remaining(2)];
% % end
% % 
% % function edge_weights = assign_tree_edge_weights_correct(n, tree)
% %     % 初始化每个节点的值为它的编号
% %     node_values = 1:n;
% %     
% %     % 初始化每条边的权值
% %     edge_weights = zeros(size(tree, 1), 1);
% %     
% %     % 逐层处理叶子节点，递归赋值边的权值
% %     while size(tree, 1) > 1
% %         [tree, node_values, edge_weights] = process_leaves_correct(tree, node_values, edge_weights);
% %     end
% %     
% %     % 如果还剩一条边，处理最后一条边
% %     if size(tree, 1) == 1
% %         remaining_edge = tree(1, :);
% %         edge_weights(end) = max(node_values(remaining_edge)); % 取两侧中较大的值
% %     end
% % end
% % 
% % function [updated_tree, updated_node_values, updated_edge_weights] = process_leaves_correct(tree, node_values, edge_weights)
% %     % 找到度为 1 的节点（叶子节点）
% %     degree = zeros(max(max(tree(:, 1:2))), 1);  % 初始化所有节点的度数
% %     for i = 1:size(tree, 1)
% %         degree(tree(i, 1)) = degree(tree(i, 1)) + 1;
% %         degree(tree(i, 2)) = degree(tree(i, 2)) + 1;
% %     end
% %     leaf_nodes = find(degree == 1);  % 度为 1 的节点即叶子节点
% %     
% %     % 处理每个叶子节点
% %     for i = 1:length(leaf_nodes)
% %         leaf = leaf_nodes(i);
% %         
% %         % 找到这个叶子节点连接的边
% %         for j = 1:size(tree, 1)
% %             if tree(j, 1) == leaf || tree(j, 2) == leaf
% %                 % 获取连接的另一个节点
% %                 connected_node = tree(j, 1) + tree(j, 2) - leaf;
% %                 
% %                 % 如果这条边没有赋过权值，给它赋值
% %                 if edge_weights(j) == 0
% %                     edge_weights(j) = node_values(leaf);  % 当前叶子的值赋给该边
% %                 end
% %                 
% %                 % 更新节点值：将叶子节点值与相连节点的值累加，更新相连节点
% %                 node_values(connected_node) = node_values(connected_node) + node_values(leaf);
% %                 
% %                 % 删除这个叶子节点和边
% %                 tree(j, :) = [];
% %                 break;
% %             end
% %         end
% %     end
% %     
% %     updated_tree = tree;
% %     updated_node_values = node_values;
% %     updated_edge_weights = edge_weights;
% % end
% % 
% % 
% % 
% %% 生成所有普吕弗序列，并将其转换为树
% function process_all_trees(n)
%     prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
%     generate_all_trees(prufer_sequence, 1, n);
% end
% %% 递归生成普吕弗序列并处理每棵树
% function generate_all_trees(seq, idx, n)
%     if idx > length(seq)
%         tree = prufer_to_tree(seq, n);
%         disp('Generated Tree:');
%         disp(tree);  % 显示生成的树（边列表）
% 
%         % 调用计算权值的函数
%         edge_weights = assign_tree_edge_weights_correct(n, tree);
%         
%         % 显示树的边及权重
%         disp('Edges with Weights:');
%         for i = 1:size(tree, 1)
%             fprintf('%5d %5d %5d\n', tree(i, 1), tree(i, 2), edge_weights(i));
%         end
%         
%         % 计算并显示每棵树的边权和
%         edge_sum = sum(edge_weights);  % 计算边权和
%         disp(['Sum of edge weights: ', num2str(edge_sum)]);
%         disp('----------------------------');
%         
%     else
%         for i = 1:n
%             seq(idx) = i;
%             generate_all_trees(seq, idx + 1, n);
%         end
%     end
% end
% %% 将普吕弗序列转换为生成树，返回树的边
% function tree = prufer_to_tree(prufer, n)
%     
%     degree = ones(1, n);  % 初始化所有节点的度数为1
%     for i = 1:length(prufer)
%         degree(prufer(i)) = degree(prufer(i)) + 1;  % 序列中的每个节点度数加1
%     end
% 
%     tree = [];  % 用于存储树的边
%     for i = 1:length(prufer)
%         for j = 1:n
%             if degree(j) == 1
%                 tree = [tree; j, prufer(i)];  % 将节点 j 和序列中的节点连接
%                 degree(j) = degree(j) - 1;
%                 degree(prufer(i)) = degree(prufer(i)) - 1;
%                 break;
%             end
%         end
%     end
% 
%     % 最后剩下两个节点，直接相连
%     remaining = find(degree == 1);
%     tree = [tree; remaining(1), remaining(2)];
% end
% 
% %% 增加映射
% function edge_weights = assign_tree_edge_weights_correct(n, tree)
%     % 初始化每个节点的值为它的编号
%     node_values = 1:n; 
%     
%     % 初始化每条边的权值
%     edge_weights = zeros(size(tree, 1), 1);
%     
%     % 复制一份树的边关系
%     remaining_tree = tree;  
%     
%     % 保存每条边在原始树中的索引
%     edge_index_map = 1:size(tree, 1); % 初始化索引映射
%     
%     % 逐层处理叶子节点
%     while size(remaining_tree, 1) > 1
%         % 计算每个节点的度数
%         degree = zeros(n, 1);
%         for i = 1:size(remaining_tree, 1)
%             degree(remaining_tree(i, 1)) = degree(remaining_tree(i, 1)) + 1;
%             degree(remaining_tree(i, 2)) = degree(remaining_tree(i, 2)) + 1;
%         end
%         
%         % 找到度为1的叶子节点
%         leaf_nodes = find(degree == 1);
%         
%         % 处理每个叶子节点
%         for i = 1:length(leaf_nodes)
%             leaf = leaf_nodes(i);
%             % 找到叶子节点连接的边
%             for j = 1:size(remaining_tree, 1)
%                 if remaining_tree(j, 1) == leaf || remaining_tree(j, 2) == leaf
%                     % 获取连接的另一个节点
%                     connected_node = remaining_tree(j, 1) + remaining_tree(j, 2) - leaf;
%                     
%                     % 给原始树中的对应边赋值，使用 edge_index_map 来找到原始索引
%                     original_edge_index = edge_index_map(j); % 通过映射找到原始索引
%                     edge_weights(original_edge_index) = node_values(leaf);
%                     
%                     % 将叶子节点的值加到相连节点上
%                     node_values(connected_node) = node_values(connected_node) + node_values(leaf);
%                     
%                     % 标记边为已处理（用NaN更清晰）
%                     remaining_tree(j, :) = NaN; 
%                     
%                     % 移除该边的映射
%                     edge_index_map(j) = NaN;
%                     break;
%                 end
%             end
%         end
%         
%         % 删除已处理的叶子节点的边
%         valid_indices = ~any(isnan(remaining_tree), 2);
%         remaining_tree = remaining_tree(valid_indices, :);
%         edge_index_map = edge_index_map(valid_indices); % 更新映射
%     end
%     
%     % 处理剩余的最后一条边
%     if size(remaining_tree, 1) == 1
%         remaining_edge = remaining_tree(1, :);
%         original_edge_index = edge_index_map(1); % 通过映射找到原始索引
%         edge_weights(original_edge_index) = max(node_values(remaining_edge)); % 取两侧节点中较大的值
%     end
% end
% 
% end

    function process_all_trees(n,x_cluster,y_cluster)
        prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
        
        % 初始化数组来存储所有树的权重和结构
        tree_weights = [];
        tree_structures = {};
        
        % 生成所有的树并存储它们的信息
        [tree_weights, tree_structures] = generate_all_trees(prufer_sequence, 1, n, tree_weights, tree_structures,x_cluster,y_cluster);
        
        % 找到最小权重的树
        [min_weight_sum, min_index] = min(tree_weights);
        T_min=min_index;
        min_cost=min_weight_sum;
        TREE=tree_structures{T_min};
        % 输出最小权重的树的结构及其权重
%         disp('Minimum Weight Tree:');
%         disp(tree_structures{T_min});  % 显示最小权值对应的树
%         disp(['Minimum Weight Sum: ', num2str(min_cost)]);  % 显示最小权值和
 % 将最小权重的树转换为0-1矩阵
    adjacency_matrix = tree_to_adjacency_matrix(TREE, n);
%  
%     % 输出0-1矩阵
%     disp('Adjacency Matrix:');
%     disp(adjacency_matrix);
    end
    
    %% 递归生成普吕弗序列并处理每棵树
    function [tree_weights, tree_structures] = generate_all_trees(seq, idx, n, tree_weights, tree_structures,x_cluster,y_cluster)
        if idx > length(seq)
            tree = prufer_to_tree(seq, n);  % 生成树
%             disp('Generated Tree:');
%             disp(tree);  % 显示生成的树（边列表）
    
            % 调用计算权值的函数
            edge_weights = assign_tree_edge_weights_correct(n, tree,x_cluster,y_cluster);
            
            % 显示树的边及权重
%             disp('Edges with Weights:');
            for i = 1:size(tree, 1)
                fprintf('%5d %5d %5d\n', tree(i, 1), tree(i, 2), edge_weights(i));
            end
            
            % 计算当前树的边权和
            edge_sum = sum(edge_weights);  % 计算边权和
%             disp(['Sum of edge weights: ', num2str(edge_sum)]);
%             disp('----------------------------');
            
            % 存储树的权重和结构
            tree_weights = [tree_weights; edge_sum];  % 将边权和加入权重数组
            tree_structures{end+1} = tree;  % 将树的结构加入元胞数组
            
        else
            for i = 1:n
                seq(idx) = i;
                [tree_weights, tree_structures] = generate_all_trees(seq, idx + 1, n, tree_weights, tree_structures,x_cluster,y_cluster);
            end
        end
    end
    
    %% 将普吕弗序列转换为生成树，返回树的边
    function tree = prufer_to_tree(prufer, n)
        
        degree = ones(1, n);  % 初始化所有节点的度数为1
        for i = 1:length(prufer)
            degree(prufer(i)) = degree(prufer(i)) + 1;  % 序列中的每个节点度数加1
        end
    
        tree = [];  % 用于存储树的边
        for i = 1:length(prufer)
            for j = 1:n
                if degree(j) == 1
                    tree = [tree; j, prufer(i)];  % 将节点 j 和序列中的节点连接
                    degree(j) = degree(j) - 1;
                    degree(prufer(i)) = degree(prufer(i)) - 1;
                    break;
                end
            end
        end
    
        % 最后剩下两个节点，直接相连
        remaining = find(degree == 1);
        tree = [tree; remaining(1), remaining(2)];
    end
    
    %% 增加映射
    function edge_weights = assign_tree_edge_weights_correct(n, tree,x_cluster,y_cluster)
        T=48;
        P_max=zeros(n,1);
        L_price=zeros(n,1);
        node_values=zeros(n,n,T);
        L_price_one=0;
        L_price_last=0;
        M1=150; % 按照最大功率划分价格等级
        M2=50;
        M3=10;
        % 初始化每个节点的值为它的功率交互矩阵
%         for i=1:n
%            node_values(i,:,:) =P(i,:,:); 
%         end
        node_values=P; 
        % 初始化建筑节点之间的距离矩阵
        distance_matrix=zeros(n,n);
        % 计算两两之间的距离
        for i = 1:n
            for j = 1:n
                distance_matrix(i, j) = sqrt((x_cluster(i) - x_cluster(j))^2 + (y_cluster(i) - y_cluster(j))^2);
            end
        end
        % 初始化每条边的权值
        edge_weights = zeros(size(tree, 1), 1);
        
        % 复制一份树的边关系
        remaining_tree = tree;  
        
        % 保存每条边在原始树中的索引
        edge_index_map = 1:size(tree, 1); % 初始化索引映射
        if size(remaining_tree, 1) == 1
            remaining_edge = remaining_tree(1, :);
            original_edge_index = edge_index_map(1); % 通过映射找到原始索引
            P_max_onedge= max(abs(sum(node_values(remaining_edge(1,1),:,:),2)));
                switch true%导线横截面积
                    case P_max_onedge>=M1
                         L_price_one=150.48;
                    case P_max_onedge>=M2&&P_max_onedge<M1
                         L_price_one=52.11; 
                    case P_max_onedge>=M3&&P_max_onedge<M2
                         L_price_one=9.77;
                    case P_max_onedge<M3
                         L_price_one=2;
                end
            edge_weights(original_edge_index) =L_price_one*distance_matrix(remaining_edge(1,1), remaining_edge(1,2)); % 取两侧节点中较大的值
        else
        % 逐层处理叶子节点
        while size(remaining_tree, 1) > 1
            % 计算每个节点的度数
            degree = zeros(n, 1);
            for i = 1:size(remaining_tree, 1)
                degree(remaining_tree(i, 1)) = degree(remaining_tree(i, 1)) + 1;
                degree(remaining_tree(i, 2)) = degree(remaining_tree(i, 2)) + 1;
            end
            
            % 找到度为1的叶子节点
            leaf_nodes = find(degree == 1);
            
            % 处理每个叶子节点
            for i = 1:length(leaf_nodes)
                leaf = leaf_nodes(i);
                % 找到叶子节点连接的边
                for j = 1:size(remaining_tree, 1)
                    if remaining_tree(j, 1) == leaf || remaining_tree(j, 2) == leaf
                        % 获取连接的另一个节点
                        connected_node = remaining_tree(j, 1) + remaining_tree(j, 2) - leaf;
                        % 给原始树中的对应边赋值，使用 edge_index_map 来找到原始索引
                        original_edge_index = edge_index_map(j); % 通过映射找到原始索引
                        P_max(leaf)= max(abs(sum(node_values(leaf,:,:),2)));
                        switch true%导线横截面积
                           case P_max(leaf)>=M1
                                 L_price(leaf)=150.48;% 【价格规格是几毛钱/米  400/10/365 待查】
                           case P_max(leaf)>=M2&&P_max(leaf)<M1
                                 L_price(leaf)=52.11; 
                           case P_max(leaf)>=M3&&P_max(leaf)<M2
                                 L_price(leaf)=9.77; 
                           case P_max(leaf)<M3
                                 L_price(leaf)=2;
                        end
                        edge_weights(original_edge_index) =L_price(leaf)*distance_matrix(leaf, connected_node);
                        
                        % 将叶子节点的值加到相连节点上
                        node_values(connected_node,:,:) = node_values(connected_node,:,:) + node_values(leaf,:,:);
                        
                        % 标记边为已处理（用NaN更清晰）
                        remaining_tree(j, :) = NaN; 
                        
                        % 移除该边的映射
                        edge_index_map(j) = NaN;
                        break;
                    end
                end
            end
            
            % 删除已处理的叶子节点的边
            valid_indices = ~any(isnan(remaining_tree), 2);
            remaining_tree = remaining_tree(valid_indices, :);
            edge_index_map = edge_index_map(valid_indices); % 更新映射
        end
        
        % 处理剩余的最后一条边
        if size(remaining_tree, 1) == 1
            remaining_edge = remaining_tree(1, :);
            original_edge_index = edge_index_map(1); % 通过映射找到原始索引
            P_max_last= max(abs(max(sum(node_values(remaining_edge(1,1),:,:),2))),max(abs(sum(node_values(remaining_edge(1,2),:,:),2))));
                switch true%导线横截面积
                   case P_max_last>=M1
                         L_price_last=150.48;
                   case P_max_last>=M2&&P_max_last<M1
                         L_price_last=52.11; 
                   case P_max_last>=M3&&P_max_last<M2
                         L_price_last=9.77; 
                   case P_max_last<M3
                         L_price_last=2;
                end
            edge_weights(original_edge_index) =L_price_last*distance_matrix(remaining_edge(1,1), remaining_edge(1,2)); % 取两侧节点中较大的值
        end
        end
    end
%% 辅助函数：将树的边结构转换为0-1邻接矩阵
function adjacency_matrix = tree_to_adjacency_matrix(tree, n)
    % 初始化 n x n 的零矩阵
    adjacency_matrix = zeros(n, n);
    
    % 遍历树的边，将相连的节点对应的矩阵位置设置为1
    for i = 1:size(tree, 1)
        node1 = tree(i, 1);
        node2 = tree(i, 2);
        adjacency_matrix(node1, node2) = 1;
        adjacency_matrix(node2, node1) = 1;  % 对称矩阵，保证无向图
    end
end
end

