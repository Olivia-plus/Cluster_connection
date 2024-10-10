% 计算适应度和更新速度和位置的函数 【将net_load替换成了load_curve,pv_curve,便于计算柔性负荷最优调度】
% 涉及到两个函数，一个是柔性负荷的优化调度函数，传输最大功率矩阵；
% 另外一个是线路传输的成本计算，涉及到最小生成树算法，树的权值为最大功率矩阵。目标是寻找到权值最小的生成树
function [fitness,trade_power,bigMatrix] = calculate_fitness(cluster_solution,load_curve,pv_curve, electricity_price,x,y,num_buildings,flexible_load_main,storage_capacity_main)
    %% 处理集群的分类结果，得到每个集群的基本情况 传入[集群划分结果，净负荷，电价，直流线路铺设成本，坐标x,y,建筑的数量，柔性负荷，储能的容量]
    %% 创建一个元胞数组来装集群的分类情况，以此来寻找对应的建筑
    % 假设 clusters 是一个行向量，每个元素表示对应建筑所属的集群编号
    % 假设 max_cluster 是集群的最大编号
    % 初始化元胞数组
    max_cluster=max(cluster_solution);
    % 初始化元胞数组
    cluster_info = cell(1, max_cluster);
    % 遍历每个集群，整理集群并放在cluster_info元胞矩阵里面
    for cluster_id = 1:max_cluster
        % 找到属于当前集群的粒子的索引,并形成一个新的数组
       indices =find(cluster_solution == cluster_id); %【如何将属于同一集群的建筑编号归到一个数组中去】
        % 将该数组添加到元胞数组中
        cluster_info{cluster_id} =indices;
    end

%% 集群归属问题 表征全部建筑的连接情况
% 初始化m x m矩阵
relationshipMatrix = zeros(num_buildings);
% 遍历每个集群，更新relationshipMatrix矩阵
for i = 1:max_cluster
    cluster = cluster_info{i};
    for j = 1:length(cluster)
        for k = j+1:length(cluster)
            relationshipMatrix(cluster(j), cluster(k)) = 1;
            relationshipMatrix(cluster(k), cluster(j)) = 1;
        end
    end
end
% 如果需要对角线为1（建筑与自身的关系），可以加上这一行
relationshipMatrix = relationshipMatrix + eye(num_buildings);
% 显示结果
% disp(relationshipMatrix);
    num_clusters =  max_cluster; % 所有的集群数量
    % 获取元胞数组中行向量的大小，获取单个集群的规模
    % sizes = cellfun(@length, cluster_info); 
    fitness_prob=zeros(num_clusters,1);% 装每次集群计算后的数值
    min_cost=zeros(num_clusters);
    PV_digest=zeros(num_clusters);
    fitness_best_matrix=cell(1,num_clusters);
    best_trade_volume_total_prob=zeros(num_clusters,1);
    fitness_connect=inf;
    %% 遍历每一个集群
    T=48;
%     net_load=load_curve-pv_curve;
    for c=1:num_clusters
        % 对进行柔性负荷调度函数传入的参数进行集群化处理
        m=length(cluster_info{c});% 矩阵的大小，建筑的数量
         if m~=0
        flexible_load=zeros(m,1);
        storage_capacity=zeros(m,1);
        load_curve_cluster=cell(1,m);
        pv_curve_cluster=cell(1,m);
        load_curve_cluster_array=zeros(m,T);
        pv_curve_cluster_array=zeros(m,T);
        x_cluster=zeros(m,1);
        y_cluster=zeros(m,1);
                for q=1:m
                    flexible_load(q)=flexible_load_main(cluster_info{c}(q));    
                    x_cluster(q)=x(cluster_info{c}(q));
                    y_cluster(q)=y(cluster_info{c}(q));
                    storage_capacity(q)=storage_capacity_main(cluster_info{c}(q));
                    load_curve_cluster{q}=load_curve{cluster_info{c}(q)};
                    pv_curve_cluster{q}=pv_curve{cluster_info{c}(q)};
                end
                 load_curve_cluster_array=cell2mat(load_curve_cluster');
                 pv_curve_cluster_array=cell2mat(pv_curve_cluster');
        % 柔性负荷调度
        [PV_digest(c),P]=FlexibleLoad(m,load_curve_cluster_array,pv_curve_cluster_array,flexible_load,storage_capacity);%【增加传出的最大交换功率矩阵,之后去掉柔性负荷和储能】
        [T,min_cost(c)]=connect_cost_min(P,m,x_cluster,y_cluster);%【x,y需要重新定义一下和处理】
%         %% 循环遍历一个集群所有可能的连接情况【这里的m可能需要修改成为对应的集群中矩阵的尺寸,已修改】
%         % 【TODO：最小生成树算法】
%        total_matrices=2^(m*(m-1)/2);% 总可能的矩阵数量
%         current_matrix=zeros(m,m);
%         best_matrix =zeros(m,m);
%         if m == 0 
%             fitness_prob(c,1) = 0; 
%             best_trade_volume_total_prob(c,1) = 0; 
%         elseif m == 1
%             fitness_prob(c,1) = inf;
%             best_trade_volume_total_prob(c,1) = 0; 
%             fitness_best_matrix{c}= best_matrix;
%         else
%             % 对于存在多种互联情况的建筑集群来说，选取成本最小的互联情况
%                 for r = 1:total_matrices
%                     % 将索引 k 转换成二进制，表示当前矩阵的状态
%                     binary = dec2bin(r, m*(m-1)/2);
%                     idx = 1;
%                     %% 生成对称矩阵
%                     for i = 1:m
%                         for j = i+1:m
%                             % 检查是否在对角线上，是则置0
%                             if i == j
%                                 current_matrix(i, j) = 0;
%                             else
%                                 % 根据二进制值设置当前位置的值
%                                 current_matrix(i, j) = str2double(binary(idx));
%                                 current_matrix(j, i) = str2double(binary(idx));
%                                 idx = idx + 1;
%                             end
%                         end
%                     end
%                     
                    % 对当前矩阵进行处理，可以在这里进行你需要的操作
                    % disp(current_matrix); % 输出所有的矩阵可能
            %       % 计算建筑之间的连接情况
            %         connectivity = zeros(num_buildings, num_buildings);% 建筑是否互联的关系矩阵
            %         for i = 1:num_buildings
            %             for j = 1:num_buildings
            %                 if cluster_solution(i) == cluster_solution(j)
            %                     connectivity(i, j) = 1; % 建筑 i 和建筑 j 属于同一集群，可以进行连接
            %                 end
            %             end
            %         end
            %    

        %% 计算电能交易量和直流线路铺设成本,这是对一个集群的，当然要遍历所有的集群【还没有遍历，已完成】这里的num——buildings需要修改为划分集群后的建筑的数量，需要重新编号
%         trade_income_total =trade_income_total+PV_digest*electricity_price;%总交易电量
%       trade_volume=zeros(1,48);%实时交易电量，临时的变量
%                     dc_cost_total=0;
%                     for i = 1:m
%                         for j = 1:m
%                             if current_matrix(i, j) == 1 % 如果两个建筑互联了【TODO：不需要这么复杂了，只需要考虑整体的建筑的交易即可，考虑跨建筑电量交易情况】
%                                 dc_cost_total=dc_cost_total+current_matrix(i, j)*dc_cost_p*sqrt((x(cluster_info{c}(i))-x(cluster_info{c}(j)))^2+(y(cluster_info{c}(i))-y(cluster_info{c}(j)))^2); % 直流线路铺设成本
%                                 for k=1:48
%                                     if (net_load{cluster_info{c}(i)}(k) < 0 && net_load{cluster_info{c}(j)}(k) >0)||(net_load{cluster_info{c}(i)}(k) > 0 && net_load{cluster_info{c}(j)}(k) <0) % 如果建筑 i,j曲线互补
%                                         trade_volume(k) = min(abs(net_load{cluster_info{c}(i)}(k)),abs(net_load{cluster_info{c}(j)}(k))); % 建筑 i 与建筑 j在k时刻交易的电量
%                                         % 更新建筑交易之后的净负荷曲线
%                                         net_load{i}(k)=net_load{cluster_info{c}(i)}(k)-sign(net_load{cluster_info{c}(i)}(k))*trade_volume(k);
%                                         net_load{j}(k)=net_load{cluster_info{c}(j)}(k)-sign(net_load{cluster_info{c}(j)}(k))*trade_volume(k);
%                                         trade_volume_total=trade_volume_total+trade_volume(k);
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                     trade_cost_total = trade_volume_total *electricity_price; % 计算交易电量收益 售电方的收益和购电方节省的电量
        % 计算适应度 目标函数
%                     total_cost =dc_cost_total - trade_cost_total; % 总成本为直流线路成本减去交易收入
%         total_cost =total_cost+min_cost; % 直流互联成本
        %加入一个判断
%                     if total_cost < fitness_connect
%                         fitness_connect = total_cost;
%                         best_matrix = current_matrix;
        best_matrix = T;
%             best_trade_volume_total=trade_volume_total;
%                     end
        fitness_best_matrix{c}= best_matrix;
%         fitness_prob(c,1) = fitness_connect; % 可能的最小总成本
         end
       
%         best_trade_volume_total_prob(c,1)= best_trade_volume_total;
%     end
% end
    end %所有集群遍历完毕，产生一种集群互联的情况，并带入下方集群模块度的计算中去
% 示例数据 归并成大的连接矩阵
% clusters = { [1, 2, 3], [4, 5] }; % 建筑编号
% connections = { [0 1 1; 1 0 1; 1 1 0], [0 1; 1 0] }; % 连接矩阵
% 计算总建筑数量
% totalBuildings = max(cellfun(@max, cluster_info));
% 计算总建筑数量
totalBuildings = 0;
cluster_info(cellfun(@isempty,cluster_info))=[];
fitness_best_matrix(cellfun(@isempty,fitness_best_matrix))=[];
for i = 1:length(cluster_info)
    totalBuildings = max(totalBuildings, max(cluster_info{i}));
end
% 初始化大矩阵
bigMatrix = zeros(totalBuildings);

% 遍历每个集群
for i = 1:length(cluster_info)
    cluster = cluster_info{i};
    connMatrix = fitness_best_matrix{i};
    % 确定集群的建筑编号
%     numBuildings = length(cluster);
    % 更新大矩阵中的连接关系
    for j = 1:length(cluster)
        for k = 1:length(cluster)
            bigMatrix(cluster(j), cluster(k)) = connMatrix(j, k);
        end
    end
end

% 确保矩阵对称
bigMatrix = (bigMatrix + bigMatrix')/2;

% e = complementarity(net_load,num_buildings);
% roh1 = modularity(relationshipMatrix,bigMatrix,e);

    %% 适应度 收益减去成本
    fitness = -(sum(PV_digest)*electricity_price-sum(min_cost));
    trade_power=sum(PV_digest);
end