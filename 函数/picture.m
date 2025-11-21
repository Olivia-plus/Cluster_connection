%% 不同储能配比下的关键指标变化
storage_ratio = [0.01, 0.05, 0.1];
net_profit = [256.4, 378.3, 473.5];
pv_local_absorption = [4854.1, 9357.8, 9357.8];
total_absorption_rate = [91.0, 100.0, 100.0];
incremental_absorption_rate = [9.7, 16.7, 16.7];
peak_load_reduction = [8.1, 11.2, 11.2];

figure;
hold on;

% 左侧纵轴
yyaxis left;
plot(storage_ratio, net_profit, '-o', 'LineWidth', 1.5, 'DisplayName', '净收益（元/天）');
hold on;
plot(storage_ratio, pv_local_absorption, '--s', 'LineWidth', 1.5, 'DisplayName', '新增光伏本地消纳量（kWh/天）');
ylabel('净收益 & 本地消纳量');
ylim([0, max(pv_local_absorption) * 1.1]);

% 右侧纵轴
yyaxis right;
plot(storage_ratio, total_absorption_rate, '-.d', 'LineWidth', 1.5, 'DisplayName', '光伏本地消总纳率(%)');
plot(storage_ratio, incremental_absorption_rate, ':^', 'LineWidth', 1.5, 'DisplayName', '新增光伏本地消纳率(%)');
plot(storage_ratio, peak_load_reduction, '-v', 'LineWidth', 1.5, 'DisplayName', '建筑尖峰负荷平均消减(%)');
ylabel('百分比 (%)');
ylim([0, max(total_absorption_rate) * 1.1]);

% 轴标签和图例
xlabel('储能配比');
title('不同储能配比下的关键指标变化');
legend('Location', 'northwest');
grid on;
hold off;

figure;

% 子图 1：净收益
subplot(3,2,1);
plot(storage_ratio, net_profit, '-o', 'LineWidth', 1.5);
title('净收益（元/天）');
xlabel('储能配比');
grid on;

% 子图 2：新增光伏本地消纳量
subplot(3,2,2);
plot(storage_ratio, pv_local_absorption, '--s', 'LineWidth', 1.5);
title('新增光伏本地消纳量（kWh/天）');
xlabel('储能配比');
grid on;

% 子图 3：光伏本地消总纳率
subplot(3,2,3);
plot(storage_ratio, total_absorption_rate, '-.d', 'LineWidth', 1.5);
title('光伏本地消总纳率(%)');
xlabel('储能配比');
grid on;

% 子图 4：新增光伏本地消纳率
subplot(3,2,4);
plot(storage_ratio, incremental_absorption_rate, ':^', 'LineWidth', 1.5);
title('新增光伏本地消纳率(%)');
xlabel('储能配比');
grid on;

% 子图 5：建筑尖峰负荷平均消减
subplot(3,2,5);
plot(storage_ratio, peak_load_reduction, '-v', 'LineWidth', 1.5);
title('建筑尖峰负荷平均消减(%)');
xlabel('储能配比');
grid on;

% 调整布局
sgtitle('不同储能配比下的关键指标变化');


%% 最优集群划分数量与收益/成本关系
cluster_num = 1:11; % 集群划分数量
benefit_cost = [0.573126383, 0.846572855, 0.979057174, 1.033302048, 1.095792736, ...
                1.15568509, 1.238666612, 1.415376296, 1.368020658, 1.263774302, 1.263409018]; % 收益/成本

% 绘制图形
figure;
plot(cluster_num, benefit_cost, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('最优集群划分数量');
ylabel('收益/成本');
title('最优集群划分数量与收益/成本关系');
% grid on;

% 显示图形
hold off;

%% 不同储能系数下的收益/成本与集群划分数量的关系
cluster_num = 4:10; % 集群划分数量
data_0 = [1.033302048, 1.095792736, 1.15568509, 1.238666612, 1.415376296, 1.368020658, 1.263774302];
data_0_01 = [1.079944043, 1.198356479, 1.3055614, 1.264086846, 1.532641186, 1.336991589, 1.635161489];
data_0_05 = [1.216079492, 1.349293468, 1.491785178, 1.54567984, 1.60132303, 1.656726433, 1.14239464];
data_0_1 = [1.305550781, 1.298685789, 1.417802771, 1.96856932, 1.73570003, 2.058798972, 1.581572616];

% 绘制图形
figure;

% 绘制四条线，分别代表不同储能系数的数据
plot(cluster_num, data_0, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', '储能系数=0');
hold on;
plot(cluster_num, data_0_01, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', '储能系数=0.01');
plot(cluster_num, data_0_05, '-^', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'DisplayName', '储能系数=0.05');
plot(cluster_num, data_0_1, '-d', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm', 'DisplayName', '储能系数=0.1');

% 添加图例、标签和标题
xlabel('集群划分数量');
ylabel('收益/成本');
title('不同储能系数下的收益/成本与集群划分数量的关系');
legend show;

% 显示网格
% grid on;

% 保持图形
hold off;

%% 储能配比与净收益的关系
storage_ratio = [0, 0.01, 0.05, 0.1]; % 储能配比
net_profit = [217.2, 283.5, 378.3, 463.0]; % 净收益（元/天）

% 绘制图形
figure;

% 绘制净收益与储能配比的关系
plot(storage_ratio, net_profit, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% 添加标签和标题
xlabel('储能配比');
ylabel('净收益（元/天）');
title('储能配比与净收益的关系');

% 显示网格
% grid on;
