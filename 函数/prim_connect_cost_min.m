function [adjacency_matrix,min_cost]=prim_connect_cost_min(P,n,x_cluster,y_cluster)
%     x = rand(nodenums, 1) * 10; % 随机生成 X 坐标（范围 0 到 10）
%     y = rand(nodenums, 1) * 10; % 随机生成 Y 坐标（范围 0 到 10）
    nodes = [x_cluster, y_cluster]; % 合并为节点坐标矩阵
    
    % 计算距离矩阵
    distMatrix = zeros(n);
    for p = 1:n
        for q= 1:n
            distMatrix(p, q) = norm(nodes(p, :) - nodes(q, :)); % 欧几里得距离
        end
    end
    
    % 使用 Prim 算法生成最小生成树
    selected = false(n, 1); % 跟踪已选择的节点
    selected(1) = true; % 从第一个节点开始
    edges = []; % 存储边
    while sum(selected) < n
        minEdge = inf; % 初始化最小边的权重
        for p = 1:n
            if selected(p) % 已选择的节点
                for q = 1:n
                    if ~selected(q) && distMatrix(p, q) < minEdge
                        minEdge = distMatrix(p, q); % 更新最小边的权重
                        u = p; % 起始节点
                        v = q; % 终止节点
                    end
                end
            end
        end
        edges = [edges; u, v]; % 存储边
        selected(v) = true; % 选择新的节点
    end
    
    % assign_weights(edges, edgeWeights, nodes);
    edge_weights=assign_weights(n,edges,distMatrix,P);
    % 计算当前树的边权和
    edge_sum = sum(edge_weights);  % 计算边权和
    min_cost = edge_sum;
   
    adjacency_matrix = tree_to_adjacency_matrix(edges, n);
    % 可视化结果
%     figure;
%     hold on;
%     for k = 1:size(edges, 1)
%         plot(nodes(edges(k, 1:2), 1), nodes(edges(k, 1:2), 2), 'b-'); % 画边
%     end
%     plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerSize', 10); % 节点
%     xlabel('X坐标');
%     ylabel('Y坐标');
%     title('随机生成节点的最短路径连接');
%     grid on;

    % 标记节点编号
    for p = 1:n
        text(nodes(p,1), nodes(p,2), num2str(p), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end

    % 显示树的边及权重
    disp('Edges with Weights:');
    for p = 1:size(edges, 1)
        fprintf('%5d %5d %5d\n', edges(p, 1), edges(p, 2), edge_weights(p));
    end
    %% 边赋值的递归函数
    function edge_weights=assign_weights(n,tree,distance_matrix,P)
        P_max=zeros(n,1);
        L_price=zeros(n,1);
        node_values=P; 
        edge_weights = zeros(size(tree, 1), 1); % 初始化边权值
        % 复制一份树的边关系
            remaining_tree = tree;  
         % 保存每条边在原始树中的索引
            edge_index_map = 1:size(tree, 1); % 初始化索引映射
            if size(remaining_tree, 1) == 1
                remaining_edge = remaining_tree(1, :);
                original_edge_index = edge_index_map(1); % 通过映射找到原始索引
                P_max_onedge= max(abs(sum(node_values(remaining_edge(1,1),:,:),2)));
                P_max_onedge=P_max_onedge*1000/(3*380*0.85);%将功率转换为孔径  
                L_price_one=P_max_onedge;
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
                            P_max(leaf)=P_max(leaf)*1000/(3*380*0.85);%将功率转换为孔径
                            L_price(leaf)=P_max(leaf);
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
                    P_max_last=P_max_last*1000/(3*380*0.85);
                    L_price_last=P_max_last;
                    edge_weights(original_edge_index) =L_price_last*distance_matrix(remaining_edge(1,1), remaining_edge(1,2)); % 取两侧节点中较大的值
                end
            end
    end

    %% 辅助函数：将树的边结构转换为0-1邻接矩阵
    function  adjacency_matrix = tree_to_adjacency_matrix(tree, n)
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