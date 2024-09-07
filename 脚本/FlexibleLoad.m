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


%% 添加24小时
%% 说明：建筑的数量是变量，时间分段暂时设置的是24小时，需要改成48小时，随机生成的实例数据需要根据自己设计的参数进行改变，储能的参数，需要根据现有的相关论文进行设计，建筑互联矩阵这个也得更改
% 参数设置
m = 5; % 建筑数量
T = 24; % 时间分段（24小时）

% 随机生成示例数据
PV_generation = rand(m, T) * 100;  % 每个建筑在每个时段的PV发电量
load_demand = rand(m, T) * 80;     % 每个建筑在每个时段的电力需求
flexible_load = rand(m, 1) * 30;   % 每个建筑的柔性负荷（可平移的负荷总量）
ev_load = rand(m, 1) * 20;         % 每个建筑的电动汽车负荷（总量）
line_capacity = rand(m, m) * 50;   % 建筑之间线路容量

% 储能系统参数
storage_capacity = rand(m, 1) * 50; % 每个建筑的储能容量
initial_soc = storage_capacity / 2; % 初始储能状态 (设置为容量的一半)
charge_rate = rand(m, 1) * 10;      % 储能充电速率
discharge_rate = rand(m, 1) * 10;   % 储能放电速率

% 互联矩阵，1表示建筑之间相互连接，0表示不连接
adjacency_matrix = [1 1 1 0 1;
                    1 1 1 0 0;
                    1 1 1 1 0;
                    0 0 1 1 1;
                    1 0 0 1 1];

% 优化变量
cvx_begin
    variable transfer(m, m, T) % 每个时段建筑之间的能量传输
    variable flex_dispatch(m, T) % 每个时段柔性负荷调度
    variable ev_dispatch(m, T) % 每个时段电动汽车负荷调度
    variable grid_feed(m, T) % 每个时段回馈给电网的能量
    variable soc(m, T) % 储能系统的状态 (State of Charge)
    variable charge(m, T) % 每个时段储能充电
    variable discharge(m, T) % 每个时段储能放电
    
    % 目标函数：最大化光伏消纳率
    maximize(sum(sum(PV_generation - grid_feed)))
    
    % 约束条件
    subject to
        % 能量平衡
        for i = 1:m
            for t = 1:T
                % 能量平衡：PV发电 + 接收能量 + 放电 = 需求 + 柔性负荷 + 电动汽车负荷 + 回馈电网的能量 + 传出能量 + 充电
                PV_generation(i, t) + sum(transfer(:, i, t)) + discharge(i, t) + flex_dispatch(i, t) + ev_dispatch(i, t) ...
                    - load_demand(i, t) == grid_feed(i, t) + sum(transfer(i, :, t)) + charge(i, t);
            end
        end
        
        % 线路容量限制
        for i = 1:m
            for j = 1:m
                if adjacency_matrix(i, j) == 0
                    transfer(i, j, :) == 0; % 没有连接的建筑间无法进行能量交易
                else
                    for t = 1:T
                        transfer(i, j, t) <= line_capacity(i, j); % 受线路容量限制
                    end
                end
            end
        end
        
        % 柔性负荷与电动汽车负荷调度限制
        for i = 1:m
            sum(flex_dispatch(i, :)) <= flexible_load(i); % 柔性负荷调度总量限制
            sum(ev_dispatch(i, :)) <= ev_load(i); % 电动汽车负荷调度总量限制
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
            soc(i, 1) == initial_soc(i); % 初始状态
            for t = 2:T
                soc(i, t) == soc(i, t-1) + charge(i, t-1) - discharge(i, t-1); % 状态平衡
            end
            soc(i, T) == initial_soc(i); % 最终状态等于初始状态
        end
        
%         % 回馈电网的能量不能为负值
%         grid_feed >= 0;

cvx_end

% 输出结果
disp('建筑之间的能量传输:')
disp(transfer)
disp('柔性负荷调度:')
disp(flex_dispatch)
disp('电动汽车负荷调度:')
disp(ev_dispatch)
disp('回馈电网的能量:')
disp(grid_feed)
disp('储能系统的状态:')
disp(soc)
disp('储能充电:')
disp(charge)
disp('储能放电:')
disp(discharge)

% 假设已经执行了优化模型，并得到了 transfer, charge, discharge, grid_feed 变量的值

% 设置颜色
colors = lines(m);

% 绘制24小时内每个建筑的能量传输情况
 figure;
for i = 1:m
    subplot(m,1,i); % 创建一个m行1列的子图布局
    bar(1:T, squeeze(sum(transfer(i,:,:), 2)), 'stacked'); % 叠加柱状图表示该建筑与其他建筑的能量传输
    hold on;
    bar(1:T, discharge(i,:), 'FaceColor', colors(i,:), 'EdgeColor', 'none'); % 叠加储能放电
    bar(1:T, -charge(i,:), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); % 叠加储能充电
    bar(1:T, grid_feed(i,:), 'FaceColor', [0 0 0], 'EdgeColor', 'none'); % 叠加与电网的能量交互
    hold off;
    title(['建筑 ' num2str(i) ' 的电能交互情况']);
    xlabel('时间 (小时)');
    ylabel('能量 (kWh)');
    legend({'建筑间传输', '储能放电', '储能充电', '与电网交互'}, 'Location', 'best');
    ylim([-max(max(max(transfer))) max(max(max(transfer)))]); % 设置y轴范围
end

% 设置整个图形的标题
sgtitle('建筑之间的能量交互情况');

