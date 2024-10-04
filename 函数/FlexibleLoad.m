% % 建筑数量
% m = 10; % 你可以根据实际建筑数量进行调整
% 
% % 参数设置
% PV_generation = rand(m, 1) * 100;  % 每个建筑的PV发电量
% load_demand = rand(m, 1) * 80;     % 每个建筑的电力需求
% flexible_load = rand(m, 1) * 30;   % 每个建筑的柔性负荷
% ev_load = rand(m, 1) * 20;         % 每个建筑的电动汽车负荷
% line_capacity = rand(m, m) * 50;   % 建筑之间线路容量
% 
% % 优化变量 CVX是一个基于matlab的凸优化建模系统
% cvx_begin
%     variable transfer(m, m) % 建筑之间的能量传输
%     variable flex_dispatch(m) % 柔性负荷调度
%     variable ev_dispatch(m) % 电动汽车负荷调度
%     
%     % 目标函数：最大化光伏消纳率
%     maximize(sum(min(PV_generation + flex_dispatch + ev_dispatch - load_demand, 0)))
%     
%     % 约束条件
%     subject to
%         % 能量平衡
%         for i = 1:m
%             PV_generation(i) + sum(transfer(:, i)) + flex_dispatch(i) + ev_dispatch(i) >= load_demand(i);
%         end
%         
%         % 线路容量限制
%         for i = 1:m
%             for j = 1:m
%                 if i ~= j
%                     transfer(i, j) <= line_capacity(i, j);
%                 end
%             end
%         end
%         
%         % 柔性负荷与电动汽车负荷调度限制
%         flex_dispatch >= 0;
%         ev_dispatch >= 0;
%         flex_dispatch <= flexible_load;
%         ev_dispatch <= ev_load;
% 
% cvx_end
% 
% % 输出结果
% disp('建筑之间的能量传输:')
% disp(transfer)
% disp('柔性负荷调度:')
% disp(flex_dispatch)
% disp('电动汽车负荷调度:')
% disp(ev_dispatch)
% 
% 
% %% CVX测试
% % clc;clear;close
% % 
% % m =20;n =10;p =4;
% % A =randn(m,n);b =randn(m,1);
% % C =randn(p,n);d =randn(p,1);e =rand;
% % cvx_begin
% %     variable x(n)
% %     minimize(norm(A *x -b,2) )
% %     subject to
% %         C*x==d;
% %         norm(x,Inf )<=e;
% % cvx_end
% 
% %% 加入储能
% % 建筑数量
% m = 5; % 你可以根据实际建筑数量进行调整
% 
% % 参数设置
% PV_generation = rand(m, 1) * 100;  % 每个建筑的PV发电量
% load_demand = rand(m, 1) * 80;     % 每个建筑的电力需求
% flexible_load = rand(m, 1) * 30;   % 每个建筑的柔性负荷
% ev_load = rand(m, 1) * 20;         % 每个建筑的电动汽车负荷
% line_capacity = rand(m, m) * 50;   % 建筑之间线路容量
% 
% % 储能系统参数
% storage_capacity = rand(m, 1) * 50; % 每个建筑的储能容量
% initial_soc = storage_capacity / 2; % 初始储能状态 (设置为容量的一半)
% charge_rate = rand(m, 1) * 10;      % 储能充电速率
% discharge_rate = rand(m, 1) * 10;   % 储能放电速率
% 
% % 互联矩阵，1表示建筑之间相互连接，0表示不连接
% adjacency_matrix = [1 1 0 0 1;
%                     1 1 1 0 0;
%                     0 1 1 1 0;
%                     0 0 1 1 1;
%                     1 0 0 1 1];
% 
% % 优化变量
% cvx_begin
%     variable transfer(m, m) % 建筑之间的能量传输
%     variable flex_dispatch(m) % 柔性负荷调度
%     variable ev_dispatch(m) % 电动汽车负荷调度
%     variable grid_feed(m) % 回馈给电网的能量
%     variable soc(m) % 储能系统的状态 (State of Charge)
%     variable charge(m) % 储能充电
%     variable discharge(m) % 储能放电
%     
%     % 目标函数：最大化光伏消纳率
%     maximize(sum(PV_generation - grid_feed))
%     
%     % 约束条件
%     subject to
%         % 能量平衡
%         for i = 1:m
%             % 能量平衡：PV发电 + 接收能量 + 放电 = 需求 + 柔性负荷 + 电动汽车负荷 + 回馈电网的能量 + 传出能量 + 充电
%             PV_generation(i) + sum(transfer(:, i)) + discharge(i) + flex_dispatch(i) + ev_dispatch(i) - load_demand(i) ...
%                 == grid_feed(i) + sum(transfer(i, :)) + charge(i);
%         end
%         
%         % 线路容量限制
%         for i = 1:m
%             for j = 1:m
%                 if adjacency_matrix(i, j) == 0
%                     transfer(i, j) == 0; % 没有连接的建筑间无法进行能量交易
%                 else
%                     transfer(i, j) <= line_capacity(i, j); % 受线路容量限制
%                 end
%             end
%         end
%         
%         % 柔性负荷与电动汽车负荷调度限制
%         flex_dispatch >= 0;
%         ev_dispatch >= 0;
%         flex_dispatch <= flexible_load;
%         ev_dispatch <= ev_load;
%         
%         % 储能充放电限制
%         charge >= 0;
%         charge <= charge_rate;
%         discharge >= 0;
%         discharge <= discharge_rate;
%         
%         % 储能容量限制
%         soc >= 0;
%         soc <= storage_capacity;
%         
%         % 储能初始和最终状态
%         soc == initial_soc + charge - discharge;
%         
%         % 回馈电网的能量不能为负值
%         grid_feed >= 0;
% 
% cvx_end
% 
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


%% 考虑柔性负荷以后，计算各集群光伏最大的消纳量，并传出，用于适应度函数的计算
%% 说明：针对一个集群而言，故建筑的数量是变量，还需要知道建筑的编号
% 储能的参数，需要根据现有的相关论文进行设计，建筑互联矩阵这个也得更改【无关互联矩阵】
% 加入线路容量 市场电缆的粗细，成本 量化方法，成本代替路径作为线路权重 分档次，土建的成本
function[y,P_transMax_array]=FlexibleLoad(Build_num,load_curve_cluster,pv_curve_cluster,flexible_load,storage_capacity)% 【TODO：不传入净负荷，改成光伏和负荷,已改】
% 参数设置 传入的是建筑的各种信息，对应编号的信息，日光伏和日净负荷曲线，可转移负荷，可平移负荷和电动汽车之类的东西，不需要互联的信息了，很好
m = Build_num; % 建筑数量【待传入】
 T = 48; % 时间分段（24小时48个点）

% 随机生成示例数据，根据建筑实际的容量配置【待传入】【可平移、可转移】
% PV_generation = rand(m, T) * 100*40;  % 每个建筑在每个时段的PV发电量
% load_demand = rand(m, T) * 80*40;     % 每个建筑在每个时段的电力需求
% flexible_load = rand(m, 1) * 30*40;   % 每个建筑的柔性负荷（可转移的负荷总量）【TODO：每个建筑设有自己的柔性负荷 1/3】
 ev_load = rand(m, 1) * 20*40*0;         % 每个建筑的电动汽车负荷（总量）
% line_capacity =ones(m,m)* 0.75*1000*135;   % 建筑之间线路容量101.25KW【由老师发的资料数据所得】
 line_capacity = ones(m, m) * 101.25; % 线路容量，每条线还不一样？怪异【暂时不考虑】
P_transMax_array = zeros(m,m);
% https://www.bilibili.com/read/cv33635168/
% 储能系统参数
initial_soc = storage_capacity *0.2; % 初始储能状态 (设置为容量的一半)
charge_rate = ones(m, 1) * 0.125;      % 储能充电速率【TODO：是否合理？一般储能电池容量和充放电速度之间的关系。0.25C用于调峰，48点，还要额外除以2。已完成】
% discharge_rate = rand(m, 1) * 0.125;   % 储能放电速率
discharge_rate = ones(m, 1) * 0.125;   % 储能放电速率
% % 互联矩阵，1表示建筑之间相互连接，0表示不连接【不需要了】
% adjacency_matrix = [1 1 1 1 1;
%                     1 1 1 1 1;
%                     1 1 1 1 1;
%                     1 1 1 1 1;
%                     1 1 1 1 1];
 
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
%     y=sum(sum((net_load_cluster>0)- grid_feed));【所有光伏量-未被消纳的（grid_feed<0）】
    y=sum((pv_curve_cluster-load_curve_cluster)>0-grid_feed>0);% 【逻辑有误,已修正】
    maximize(y)
    
    % 约束条件
    subject to
        % 能量平衡
        for i = 1:m
            for t = 1:T
                % 能量平衡：PV发电 + 接收能量 + 放电 = 需求 + 柔性负荷 + 电动汽车负荷 + 回馈电网的能量 + 传出能量 + 充电
              pv_curve_cluster{i}(t)+ sum(transfer(:, i, t)) + discharge(i, t)*storage_capacity(m)+ flex_dispatch(i, t) + ev_dispatch(i, t) ...
                    == grid_feed(i, t) + sum(transfer(i, :, t)) + charge(i, t)*storage_capacity(m)+load_curve_cluster{i}(t);
            end
        end
        
        % 线路容量限制
        for i = 1:m
            for j = 1:m
%                 if adjacency_matrix(i, j) == 0
%                     transfer(i, j, :) == 0; % 没有连接的建筑间无法进行能量交易【现在可以进行交易了】
%                 else
                    for t = 1:T
                        transfer(i, j, t) <= line_capacity(i, j); % 受线路容量限制 其实不是很受限制才对
                    end
                    P_transMax_array(i,j)= sum(transfer(i, j, :));
%                 end
            end
        end
        
        % 柔性负荷与电动汽车负荷调度限制
        for i = 1:m
            sum(flex_dispatch(i, :)) = flexible_load(i); % 柔性负荷调度总量限制
            sum(ev_dispatch(i, :)) = ev_load(i); % 电动汽车负荷调度总量限制
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
            soc(i, :) >= 0;
            soc(i, :) <= storage_capacity(i);
        end
        
        % 储能状态平衡和初始条件
        for i = 1:m
            soc(i, 1) = initial_soc(i); % 初始状态
            for t = 2:T
                soc(i, t) = soc(i, t-1) + charge(i, t-1)*storage_capacity(i)- discharge(i, t-1)*storage_capacity(i); % 状态平衡
            end
            soc(i, T) = initial_soc(i); % 最终状态等于初始状态
        end
        
%         % 回馈电网的能量不能为负值 【可以为负，为负的部分就是从电网取电的过程，当光伏不能满足负荷要求的时候，就会向电网取电】
%         grid_feed >= 0;

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
