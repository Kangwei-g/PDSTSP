function [center,guilei]=K_mean(k,Distance_1,distance_min_1)

k ; %聚类数量


point_num=size(Distance_1,1);
point_num; % d为样本点的个数

%第一步：初始随机选择k个点作为类心。
%第一步分析：先随机找一个解，作为初始可行解。故初始解会影响聚类的速度和最终聚类结果（陷入局部最优解）。可以通过多次随机选取初始解找到较优的结果。

center = zeros(k,1);
suiji = randperm(point_num);

% center = suiji(1:k)';
center = [6,7];


ddcs = 50; %设置迭代次数，若迭代超过一定次数则停止
zjl = zeros(1,ddcs); %总距离，所有类心到类内各点距离之和。

for q = 1:1:ddcs
    %第二步，将所有点，按距离分到最近的类心里，得到k个类。
    %第二步分析：样本点归在离的最近的类心。
    
    %计算每个点到类心的距离
    distences = zeros(k,point_num);
    dist=distance_min_1;
    dist(dist==inf)=0;
    for i = 1:1:k
	for j = 1:1:point_num
		distences(i,j) = dist(center(i),j);
	end
    end
    %将到类心最短距离的点归类
    guilei = zeros(1,point_num);
    for i = 1:1:point_num
        guilei(1,i) = min(find(distences(:,i) == min(distences(:,i))));
    end
    
    
    
    
    %第三步：类内寻找新类心
    %第三步分析，对于一片数据点来说，数据点的中心（坐标的平均值），离这些点的总距离最短。

    %新类心到原来类内所有点距离最短  
    DistSum=[];
    for i = 1:k
        a = find(guilei(1,:) == i);
        for j=1:length(a)
            distsum=0;
            for m=1:length(a)
                distsum=distsum+dist(a(j),a(m));
            end
            DistSum(i,j)=distsum;
        end
        DistSum(DistSum==0)=inf;
        center(i)=a(find(DistSum(i,:)==min(DistSum(i,:)),1));
        
    end
    DistSum(DistSum==inf)=0;
    %第四步：若类心没有改变，则聚类结束
    %计算类心离类内点的总距离
    %第四步分析：若类心改变，则总距离变化，以总距离的变化来判断类心是否有变化。
    zjl(1,q) = sum(sum(DistSum));
    if q > 1
        if zjl(1,q) == zjl(1,q - 1)
            break
        end
    end
    
end