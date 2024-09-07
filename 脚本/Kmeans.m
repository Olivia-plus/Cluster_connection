clear;clc;
data=xlsread('C:\Users\WuAoli\Desktop\Kmeans.xlsx');
data=data(:,2:7);
%% 原理推导K均值
[m,n]=size(data); %读取数据的行数与列数
cluster_num=3; %自定义分类数
cluster=data(randperm(m,cluster_num),:);
epoch_max=1000;%最大次数
therad_lim=0.001;%中心变化阈值
epoch_num=0;
while(epoch_num<epoch_max)
    epoch_num=epoch_num+1;
    for i=1:cluster_num
        distance=(data-repmat(cluster(i,:),m,1)).^2;
    distance1(:,i)=sqrt(sum(distance'));
    end
    [~,index_cluster]=min(distance1');
    for j=1:cluster_num
    cluster_new(j,:)=mean(data(find(index_cluster==j),:));
    end
    if (sqrt(sum((cluster_new-cluster).^2))>therad_lim)
        cluster=cluster_new;
    else
        break;
    end
end
%% 画出聚类效果
figure
subplot(2,1,1) %画子图，在这里是一图上可画两张子图
a=unique(index_cluster); %找出分类出的个数
C=cell(1,length(a));
for i=1:length(a)
   C(1,i)={find(index_cluster==a(i))};
end
for j=1:cluster_num
    data_get=data(C{1,j},:);
    scatter(data_get(:,1),data_get(:,2),100,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.9);
    hold on
end
plot(cluster(:,1),cluster(:,2),'kd','LineWidth',2);
hold on
sc_t=mean(silhouette(data,index_cluster'));
title_str=['原理推导K均值聚类','  聚类数为：',num2str(cluster_num),'  SC轮廓系数:',num2str(sc_t)];
title(title_str)
%% MATLAB自带kmeans函数
subplot(2,1,2) %画子图，在这里是一图上可画两张子图
cluster_num=3; %自定义分类数
[index_km,center_km]=kmeans(data,cluster_num) ;%MATLAB自带kmeans函数
a=unique(index_km); %找出分类出的个数
C=cell(1,length(a));
for i=1:length(a)
   C(1,i)={find(index_km==a(i))};
end
for j=1:cluster_num
    data_get=data(C{1,j},:);
    scatter(data_get(:,1),data_get(:,2),100,'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.9);
    hold on
end
plot(center_km(:,1),center_km(:,2),'kd','LineWidth',2);
hold on
sc_k=mean(silhouette(data,index_km));
title_str1=['MATLAB自带kmeans函数','  聚类数为：',num2str(cluster_num),'  SC轮廓系数:',num2str(sc_k)];
title(title_str1)
