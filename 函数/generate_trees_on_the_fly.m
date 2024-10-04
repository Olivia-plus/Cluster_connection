%% 产生全部的生成树 输入节点的个数
function generate_trees_on_the_fly(n)
    % 如果只有2个节点，生成树只能是单一的一条边
    if n == 2
        disp([1, 2]);
        return;
    end

    % 逐步生成普吕弗序列，生成完后直接转换为树并打印
    prufer_sequence = zeros(1, n-2);  % 初始化长度为 n-2 的普吕弗序列
    generate_sequence(prufer_sequence, 1, n); % 递归
end

function generate_sequence(seq, idx, n)
    % 递归生成普吕弗序列并打印生成树
    if idx > length(seq)
        tree = prufer_to_tree(seq, n);
        disp(tree);
    else
        for i = 1:n
            seq(idx) = i;
            generate_sequence(seq, idx + 1, n);
        end
    end
end

function tree = prufer_to_tree(prufer, n)
    % 将普吕弗序列转换为生成树
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

% 每个时刻单向传输
% 3和4当成一个整体，负荷和光伏曲线相加，其实就是需要把他们的功率矩阵（带时间线的）相加起来就可以。
% 整个流程图 流程子图 配一个图一段画去描述
% 线径和成本计算方法 表

