% %% 绘制每个建筑的光伏发电量、固定负荷、柔性负荷和电动汽车负荷的叠加柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     % 生成叠加柱状图数据
%     bar_data = [PV_generation(i, :); load_demand(i, :); flexible_load(i, :); ev_load(i, :)];
%     bar(1:T, bar_data', 'stacked'); % 叠加柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i)]);
%     legend('PV Generation', 'Load Demand', 'Flexible Load', 'EV Load');
%     grid on;
% end
% 
% % 绘制建筑之间的能量传输的叠加柱状图
% figure;
% for t = 1:T
%     subplot(6, 4, t); % 24小时分成6行4列的子图
%     bar(1:m, squeeze(transfer(:, :, t))', 'stacked'); % 叠加柱状图
%     xlabel('Building');
%     ylabel('Energy Transfer (kW)');
%     title(['Energy Transfer at Hour ' num2str(t)]);
%     legend('Building 1', 'Building 2', 'Building 3', 'Building 4', 'Building 5');
%     grid on;
% end
% 
% % 绘制储能系统的状态 (SOC) 的叠加柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar(1:T, soc(i, :), 'FaceColor', 'k', 'DisplayName', 'State of Charge (SOC)');
%     xlabel('Hour');
%     ylabel('SOC (kWh)');
%     title(['Building ' num2str(i) ' State of Charge']);
%     grid on;
% end
% 
% % 绘制柔性负荷调度和电动汽车负荷调度的叠加柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar_data_dispatch = [flex_dispatch(i, :); ev_dispatch(i, :)];
%     bar(1:T, bar_data_dispatch', 'stacked'); % 叠加柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i) ' Load Dispatch']);
%     legend('Flexible Load Dispatch', 'EV Load Dispatch');
%     grid on;
% end
% 
% % 绘制回馈电网的能量和从电网购电的能量的叠加柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar_data_grid = [grid_feed(i, :); grid_purchase(i, :)];
%     bar(1:T, bar_data_grid', 'stacked'); % 叠加柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i) ' Grid Interaction']);
%     legend('Grid Feed', 'Grid Purchase');
%     grid on;
% end
% 
% %% 绘制每个建筑的光伏发电量、固定负荷、柔性负荷和电动汽车负荷的分离柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     % 生成分离柱状图数据
%     bar_data = [PV_generation(i, :); load_demand(i, :); flexible_load(i, :); ev_load(i, :)];
%     bar(1:T, bar_data', 'grouped'); % 分离柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i)]);
%     legend('PV Generation', 'Load Demand', 'Flexible Load', 'EV Load');
%     grid on;
% end
% 
% % 绘制建筑之间的能量传输的分离柱状图
% figure;
% for t = 1:T
%     subplot(6, 4, t); % 24小时分成6行4列的子图
%     bar(1:m, squeeze(transfer(:, :, t))', 'grouped'); % 分离柱状图
%     xlabel('Building');
%     ylabel('Energy Transfer (kW)');
%     title(['Energy Transfer at Hour ' num2str(t)]);
%     legend('To Building 1', 'To Building 2', 'To Building 3', 'To Building 4', 'To Building 5');
%     grid on;
% end
% 
% % 绘制储能系统的状态 (SOC) 的分离柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar(1:T, soc(i, :), 'FaceColor', 'k', 'DisplayName', 'State of Charge (SOC)');
%     xlabel('Hour');
%     ylabel('SOC (kWh)');
%     title(['Building ' num2str(i) ' State of Charge']);
%     grid on;
% end
% 
% % 绘制柔性负荷调度和电动汽车负荷调度的分离柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar_data_dispatch = [flex_dispatch(i, :); ev_dispatch(i, :)];
%     bar(1:T, bar_data_dispatch', 'grouped'); % 分离柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i) ' Load Dispatch']);
%     legend('Flexible Load Dispatch', 'EV Load Dispatch');
%     grid on;
% end
% 
% % 绘制回馈电网的能量和从电网购电的能量的分离柱状图
% figure;
% for i = 1:m
%     subplot(m, 1, i);
%     bar_data_grid = [grid_feed(i, :); grid_purchase(i, :)];
%     bar(1:T, bar_data_grid', 'grouped'); % 分离柱状图
%     xlabel('Hour');
%     ylabel('Power (kW)');
%     title(['Building ' num2str(i) ' Grid Interaction']);
%     legend('Grid Feed', 'Grid Purchase');
%     grid on;
% end

%% 考虑柔性负荷以后，计算各集群光伏最大的消纳量，并传出，用于适应度函数的计算
%% 说明：针对一个集群而言，故建筑的数量是变量，还需要知道建筑的编号
% 储能的参数，需要根据现有的相关论文进行设计，建筑互联矩阵这个也得更改【无关互联矩阵】
% 加入线路容量 市场电缆的粗细，成本 量化方法，成本代替路径作为线路权重 分档次，土建的成本
function[y,P_transMax_array]=FlexibleLoad(Build_num,load_curve_cluster,pv_curve_cluster,flexible_load,storage_capacity)
% 参数设置 传入的是建筑的各种信息，对应编号的信息，日光伏和日净负荷曲线，可转移负荷，可平移负荷和电动汽车之类的东西，不需要互联的信息了，很好
m = Build_num; % 建筑数量【待传入】
T = 48; % 时间分段（24小时48个点）

% 随机生成示例数据，根据建筑实际的容量配置【待传入】【可平移、可转移】
ev_load = rand(m, 1) * 20*40*0;         % 每个建筑的电动汽车负荷（总量）

% 储能系统参数
initial_soc = storage_capacity *0.5; % 初始储能状态 (设置为容量的一半)
charge_rate = ones(m, 1) * 0.125;      % 储能充电速率【TODO：是否合理？一般储能电池容量和充放电速度之间的关系。0.25C用于调峰，48点，还要额外除以2。已完成】
discharge_rate = ones(m, 1) * 0.125;   % 储能放电速率

% 优化变量
cvx_begin
    variable transfer(m, m, T) % 每个时段建筑之间的能量传输
    variable flex_dispatch(m, T) % 每个时段柔性负荷调度
    variable ev_dispatch(m, T) % 每个时段电动汽车负荷调度
    variable grid_feed(m, T) % 每个时段电网的能量【有正有负】
    variable soc(m, T) % 储能系统的状态 (State of Charge)
    variable charge(m, T) % 每个时段储能充电
    variable discharge(m, T) % 每个时段储能放电
    
    % 目标函数：最大化光伏消纳率

%%计算不带柔性负载和储能的光伏消纳量    
dragon=pv_curve_cluster- load_curve_cluster;

% 光伏总量sum(pv_curve_cluster)  
% 光伏建筑剩余光伏消纳量sum(max(pv_curve_cluster - load_curve_cluster, 0))  
% 光伏新增消纳量sum(sum(max(pv_curve_cluster -load_curve_cluster, 0)- max(grid_feed, 0))) 
% 光伏自消纳率1-sum(max(pv_curve_cluster - load_curve_cluster, 0)) /sum(pv_curve_cluster)  
% 光伏额外消纳率sum(sum(max(pv_curve_cluster -load_curve_cluster, 0)- max(grid_feed, 0))) /sum(pv_curve_cluster)
% 光伏总消纳率1-sum(max(pv_curve_cluster - load_curve_cluster, 0)) /sum(pv_curve_cluster)+sum(sum(max(pv_curve_cluster -load_curve_cluster, 0)- max(grid_feed, 0))) /sum(pv_curve_cluster) 
y =sum(sum(max(pv_curve_cluster - load_curve_cluster, 0)- max(grid_feed, 0)));
maximize(y)
    % 约束条件
    subject to
        % 能量平衡
        for i = 1:m
            for t = 1:T
                % 能量平衡：PV发电 + 接收能量 + 放电 = 馈电需求 + 柔性负荷 + 电动汽车负荷 + 回馈电网的能量 + 传出能量 + 充电
            pv_curve_cluster(i,t)+ sum(transfer(:, i, t)) + discharge(i, t)*storage_capacity(i)  ...
                    == grid_feed(i, t) + sum(transfer(i, :, t)) + charge(i, t)*storage_capacity(i)+load_curve_cluster(i,t)+ ev_dispatch(i, t)+ flex_dispatch(i, t);
            end
        end
        
       % grid_feed需要约束，即所有建筑缺电时不向电网输电，所有建筑不缺电时grid_feed不向电网索电  
       for i=1:T
           if all(dragon(:,i)<0)
               grid_feed(:,i)<=0;
           end
            if all(dragon(:,i)>0)
               grid_feed(:,i)>=0;
           end
       end

      for k=1:T
        for i=1:m
            for j=1:m
                transfer(i, j, k)+transfer(j, i, k)==0;
            end
        end
      end
        
        % 柔性负荷与电动汽车负荷调度限制
        for i = 1:m
             sum(abs(flex_dispatch(i, :))) <= flexible_load(i); % 柔性负荷调度总量限制
             sum(ev_dispatch(i, :)) == ev_load(i); % 电动汽车负荷调度总量限制【小于等于还是等于】
             flex_dispatch(i, :) >= 0;
             ev_dispatch(i, :) >= 0;
        end

          % 储能充放电限制
        for i = 1:m
            charge(i, :) >= 0;
            charge(i, :) <= charge_rate(i);
            discharge(i, :) >= 0;
            discharge(i, :) <= discharge_rate(i);
        end
    
    % 储能容量限制
    for i = 1:m
            soc(i, :) >= 0.1*storage_capacity(i);
            soc(i, :) <= 0.9*storage_capacity(i);
    end

    % 储能初始状态与状态更新（每小时）
    soc(:, 1) == initial_soc; % 初始时刻的储能状态
    for t = 2:T
        soc(:, t) == soc(:, t-1) + charge(:, t-1)*storage_capacity(i) - discharge(:, t-1)*storage_capacity(i); % 每小时的SOC变化
    end
    soc(:, T) == initial_soc; % 最终状态等于初始状态
    
cvx_end

% % 输出结果
% disp('建筑之间的能量传输:')
% disp(transfer)
% disp('柔性负荷调度:')
% disp(flex_dispatch)
% disp('电动汽车负荷调度:')
% disp(ev_dispatch)
% disp('回馈电网的能量:')
% disp(grid_feed)
% disp('储能系统的状态:')
% disp(soc)
% disp('储能充电:')
% disp(charge)
% disp('储能放电:')
% disp(discharge)

% disp('Optimized Grid Feed:');
% disp(grid_feed);
% disp('Optimized Transfer:');
% disp(transfer);

% 尖峰负荷削减分为两部分：第一部分，建筑直接互联功率互济产生的负荷削减；第二部分，加入储能之后，对负荷的进一步削减，
% load_cut_curve=load_curve_cluster-

% 计算有储能和光伏的负荷削减
% load_curve_cluster 表示建筑的初始负荷需求
% grid_feed 表示优化后从电网获取的能量
% storage_discharge 表示储能放电量
% % 1. 计算储能优化后的净负荷曲线
net_load_no_storage = max(load_curve_cluster - pv_curve_cluster, 0); % 无储能情况下的正净负荷
net_load_with_storage = max(load_curve_cluster - pv_curve_cluster, 0) - discharge;
net_load_with_storage(net_load_with_storage < 0) = 0; % 负值修正为0
% % 2. 设置绘图
T = 48; % 每天48个时间点（每30分钟一个点）
time = 1:T; % 时间向量
% 计算无储能和光伏的负荷需求，即 load_curve_cluster 原始负荷需求
load_without_pv_storage = load_curve_cluster; % 没有光伏和储能时的负荷需求

% 计算有储能和光伏的负荷需求，即考虑光伏和储能优化后的负荷
load_with_pv_storage = load_curve_cluster - (pv_curve_cluster + discharge .* storage_capacity); % 有储能和光伏后的负荷

% 负荷削减量
% load_reduction = max(0, load_without_pv_storage - load_with_pv_storage); % 防止出现负值
load_reduction = max(0, net_load_no_storage - load_with_pv_storage); % 加入互联储能之后的负荷削减


% 绘制时序柱状图
figure;
hold on;
for i = 1:m
    % 每个建筑的负荷需求柱状图，按时间分布
    subplot(m, 1, i);
    
    % 绘制无储能和光伏时的负荷
%     bar(1:T, load_without_pv_storage(i, :), 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none'); 
    bar(1:T, net_load_no_storage(i, :), 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none'); 
    hold on;
    
    % 绘制优化后的负荷需求
    bar(1:T, load_with_pv_storage(i, :), 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none');
    % 绘制储能优化后的净负荷曲线
%     plot(time, net_load_with_storage, '-', 'Color', [0, 0.45, 0.74], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '储能优化后净负荷');
    
    % 绘制负荷削减
    bar(1:T, load_reduction(i, :), 'FaceColor', [0.9290, 0.6940, 0.1250], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
    
    title(['建筑 ' num2str(i) ' 的负荷削减时序图']);
    xlabel('时间 (半小时)');
    ylabel('负荷需求 (kWh)');
    legend({'净负荷', '有储能优化', '引入储能后的负荷削减'}, 'Location', 'best');
end

% 设置整个图形的标题
sgtitle('所有建筑在各时间段的负荷削减时序图');
hold off;


% 计算每个时段的削减百分比
reduction_percentage = zeros(m, T); % 初始化削减百分比矩阵
for i = 1:m
    for t = 1:T
        if net_load_no_storage(i, t) > 0
            reduction_percentage(i, t) = (net_load_no_storage(i, t) - net_load_with_storage(i, t)) / net_load_no_storage(i, t) * 100;
        end
    end
end

% 输出结果
disp('无储能情况下的正净负荷：');
disp(net_load_no_storage);
disp('储能优化后的正净负荷：');
disp(net_load_with_storage);
disp('每个时段的负荷削减百分比：');
disp(reduction_percentage);

% Step 2: 计算每个建筑在峰值时段的削减百分比
peak_reduction_percentage = zeros(m, 1); % 初始化峰值削减百分比向量
for i = 1:m
    % 找到每个建筑的净负荷峰值及其对应的时段
    [peak_load, peak_time] = max(net_load_no_storage(i, :));
    
    % 计算储能优化后在峰值时段的削减百分比
    if peak_load > 0
        peak_reduction_percentage(i) = (peak_load - net_load_with_storage(i, peak_time)) / peak_load * 100;
    else
        peak_reduction_percentage(i) = 0; % 无峰值时，削减率设为0
    end
end

% 输出结果
disp('每个建筑的净负荷峰值削减百分比：');
disp(peak_reduction_percentage);


% 示例数据（假设已计算并导入的数据）
% net_load_no_storage: 无储能优化时的净负荷（负荷 - 光伏 > 0 的部分）
% discharge: 储能系统在每个时间段的放电量
% grid_feed: 电网交互（作为柱状图的其中一部分）
% 其他可视化数据（如：transfer等）
% 
% % 1. 计算储能优化后的净负荷曲线
% net_load_with_storage = max(load_curve_cluster - pv_curve_cluster, 0) - discharge;
% net_load_with_storage(net_load_with_storage < 0) = 0; % 负值修正为0
% 
% % 2. 设置绘图
% T = 48; % 每天48个时间点（每30分钟一个点）
% time = 1:T; % 时间向量
% 
% % 3. 绘制柱状图与曲线
% figure;
% hold on;
% 
% % 绘制无储能优化的净负荷柱状图
% bar(time, max(load_curve_cluster - pv_curve_cluster, 0), 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none', 'DisplayName', '无储能净负荷');
% 
% % 绘制储能优化后的净负荷曲线
% plot(time + 0.5, net_load_with_storage, '-o', 'Color', [0, 0.45, 0.74], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '储能优化后净负荷');
% 
% % 添加电网交互（柱状图，设置不同颜色以区分）
% bar(time, grid_feed, 'FaceColor', [0.3, 0.3, 0.3], 'EdgeColor', 'none', 'DisplayName', '电网交互');
% 
% % 添加图例和标签
% legend();
% xlabel('时间 (30分钟为单位)');
% ylabel('净负荷 (kWh)');
% title('储能优化后的负荷削减效果');
% grid on;
% 
% hold off;






% 
% % 设置经典的三种颜色
% color_storage_discharge = [0, 0.4470, 0.7410]; % 蓝色
% color_storage_charge = [0.8500, 0.3250, 0.0980]; % 橙色
% color_grid_feed = [0.9290, 0.6940, 0.1250]; % 黄色
% 
% % 绘制24小时内每个建筑的能量交互情况
% figure;
% for i = 1:m
%     subplot(m,1,i); % 创建一个m行1列的子图布局
%     bar(1:T, squeeze(sum(transfer(i,:,:), 2)), 'FaceColor', [0.3, 0.3, 0.3], 'EdgeColor', 'none'); % 灰色表示建筑间传输
%     hold on;
%     bar(1:T, discharge(i,:), 'FaceColor', color_storage_discharge, 'EdgeColor', 'none'); % 蓝色代表储能放电
%     bar(1:T, -charge(i,:), 'FaceColor', color_storage_charge, 'EdgeColor', 'none'); % 橙色代表储能充电
%     bar(1:T, grid_feed(i,:), 'FaceColor', color_grid_feed, 'EdgeColor', 'none'); % 黄色代表与电网的能量交互
%     hold off;
%     title(['建筑 ' num2str(i) ' 的电能交互情况']);
%     xlabel('时间 (小时)');
%     ylabel('能量 (kWh)');
%     legend({'建筑间传输', '储能放电', '储能充电', '与电网交互'}, 'Location', 'best');
% end
% 
% % 设置整个图形的标题
% sgtitle('建筑之间的能量交互情况');

% % 无储能和光伏情况下的负荷需求
% load_without_pv_storage = sum(load_curve_cluster, 2); % 没有光伏和储能时的总负荷需求
% 
% % 有储能和光伏优化后的负荷需求
% load_with_pv_storage = sum(load_curve_cluster - (pv_curve_cluster + discharge .* storage_capacity), 2); % 考虑光伏和储能后的负荷需求
% 
% % 负荷削减量
% load_reduction = load_without_pv_storage - load_with_pv_storage;
% 
% % 绘制柱状图
% figure;
% bar([load_without_pv_storage, load_with_pv_storage, load_reduction], 'stacked'); 
% title('储能和光伏对负荷削减的影响');
% xlabel('建筑编号');
% ylabel('负荷需求 (kWh)');
% legend({'无储能和光伏', '有储能和光伏优化', '负荷削减'}, 'Location', 'best');
% grid on;

P_transMax_array=transfer;
% 假设已经执行了优化模型，并得到了 transfer, charge, discharge, grid_feed 变量的值

% 设置颜色
% colors = lines(m);

% % 绘制24小时内每个建筑的能量传输情况
%  figure;
% for i = 1:m
%     subplot(m,1,i); % 创建一个m行1列的子图布局
%     bar(1:T, squeeze(sum(transfer(i,:,:), 2)), 'stacked'); % 叠加柱状图表示该建筑与其他建筑的能量传输
%     hold on;
%     bar(1:T, discharge(i,:), 'FaceColor', colors(i,:), 'EdgeColor', 'none'); % 叠加储能放电
%     bar(1:T, -charge(i,:), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); % 叠加储能充电
%     bar(1:T, grid_feed(i,:), 'FaceColor', [0 0 0], 'EdgeColor', 'none'); % 叠加与电网的能量交互
%     hold off;
%     title(['建筑 ' num2str(i) ' 的电能交互情况']);
%     xlabel('时间 (小时)');
%     ylabel('能量 (kWh)');
%     legend({'建筑间传输', '储能放电', '储能充电', '与电网交互'}, 'Location', 'best');
% %     ylim([-max(max(max(transfer))) max(max(max(transfer)))]); % 设置y轴范围
% end
% 
% % 设置整个图形的标题
% sgtitle('建筑之间的能量交互情况');
end

