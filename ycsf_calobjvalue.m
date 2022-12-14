% 2.2.3 计算目标函数值
% calobjvalue.m函数的功能是实现目标函数的计算，其公式采用本文示例仿真，可根据不同优化问题予以修改。
%遗传算法子程序
%Name: calobjvalue.m
%实现目标函数的计算，计算最短的时间。


function [objvalue,dist_vehicle,path_vehicle,path_UAV]=ycsf_calobjvalue(pop)
% 加载运行函数所需数据
load('data_of_2022527B_ques2.mat')
load('ques2_calculated_data.mat')

objvalue=zeros(size(pop,1),1); % 即最短配送时间
dist_vehicle=zeros(size(pop,1),1); % 车辆最优配送方案的总里程
path_vehicle=cell(size(pop,1),1); % 即车辆最优配送方案
path_UAV=cell(size(pop,1),1); % 即无人机最优配送方案

for i = 1:size(pop,1)  % 遍历所有种群，计算最短时间
    
    %% 先计算需要车辆配送的地点 VehicleDelivery 和需要无人机配送的地点 UAVDelivery
    
    a=UAVDeliveryPositionMaybe(pop(i,:)==0);
    
    b=UAVNoFlyPosition;
    
    % 需要车辆配送的地点 VehicleDelivery
    VehicleDelivery=sort([a,b]);
    % 需要无人机配送的地点 UAVDelivery
    UAVDelivery=UAVDeliveryPositionMaybe(pop(i,:)==1);
    
    %% 再使用与第一问相同的遗传算法（mtspf_ga.m函数）计算车辆的最佳配送路线
    n = length(VehicleDelivery); % 车辆配送地点数量
    xy = rand(n,2); % 调用函数这个变量用不到，因为我们已经有了两两之间距离，但函数运行仍需要输入
    salesmen = 1;  % 旅行者，即配送车辆数
    min_tour = n; % 一个车辆最少可以配送几个地点
    pop_size = 80; % 种群数
    num_iter = 5e1;  % 迭代数
    
    dmat=distance_min_1(VehicleDelivery,VehicleDelivery);
    [opt_rte,~,~] = mtspf_ga(xy,dmat,salesmen,min_tour,pop_size,num_iter,0,0);
    % opt_rte，不包括起止点，且起点和终点为距离矩阵中第一个点。
    % 还原路径原先的序号
    opt_rte_restore=VehicleDelivery(opt_rte);
    opt_rte=opt_rte_restore;
    
    % 对路径进行完善补充
    path=[VehicleDelivery(1),opt_rte]; % 将起点1补充入内
    % 又因题目要求从物资集中地出发回到物资集中地，且以上路径本身属于循环，可以从任意地点开始且回到该地点，所以使路径从物资集中地开始
    a=find(path==DistributionCenter,1);
    path=[path(a:end),path(1:a-1),DistributionCenter];
    % 在这个路径中不会出现重复的点，但是实际上有的两个点之间没有直接连通，需要先到其他点，例如6-4-6和5-2-3-5这两个回路
    % 所以对于没有直接通路的路径需要补充间接通路
    flag=1;
    while flag
        t=0;
        for j=1:length(path)-1
            if Distance_1(path(j),path(j+1))==inf
                path=[path(1:j-1),path_min_1{path(j),path(j+1)},path(j+2:end)];
                t=t+1;
            end
        end
        if t==0
            flag=0;
        end
    end
    % path=[9,8,7,5,2,5,6,10,9];
    % 记录车辆最优配送方案
    path_vehicle{i}=path;
    
    % 最佳配送路线的总里程
    distance_total=0;
    for j=1:length(path)-1
        distance_total=distance_total+distance_min_1(path(j),path(j+1));
    end
    % 记录车辆最优配送方案的总里程
    dist_vehicle(i)=distance_total;
    
    % 最佳配送路线所花时间
    time_total=distance_total/VehicleSpeed;
    
    %% 在车辆的最佳配送路线基础上分配无人机配送任务
    % 需要无人机配送的地点 UAVDelivery
    record=[];
    for j=1:length(UAVDelivery)
        % 如果无人机配送的地点在最佳配送路线中出现，则此地点无需无人机配送，修改种群中的个体
        if ~isempty(find(path==UAVDelivery(j), 1))
            record=[record,j];
        end
    end
    UAVDelivery(record)=[];
    
    % 一些地点无需无人机配送，修改该个体的染色体
    pop_change=pop(i,:);
    for j=1:length(UAVDeliveryPositionMaybe)
        if ~isempty(find(UAVDeliveryPositionMaybe(j)==UAVDelivery, 1))
            pop_change(j)=1;
        else
            pop_change(j)=0;
        end
    end
    pop(i,:)=pop_change;
    
    % 初始化包含本次无人机配送的地点的无人机路径
    path_mid=[UAVDelivery',zeros(length(UAVDelivery),size(path_end,2))];
    % 初始化无人机各配送地点的最佳路径
    path_UAV_mid=[UAVDelivery',zeros(length(UAVDelivery),2)];
    % 初始化无人机各配送地点的最佳路径距离
    dist_UAV_mid=[UAVDelivery',zeros(length(UAVDelivery),1)];
    % 初始化无人机各配送最佳路径耗时
    time_UAV_mid=[UAVDelivery',zeros(length(UAVDelivery),1)];
    
    % 对物资需求合并后仍在无人机最大载重量以内，且距离相邻的无人机进行合并
    combine=[];
    for j=1:length(UAVDelivery)
        for k=j:length(UAVDelivery)
            % 判断物资需求能否满足 且满足相邻条件
            if MaterialNeed(UAVDelivery(j))+MaterialNeed(UAVDelivery(k))<=UAVCapacity ...
                    && Distance_2(UAVDelivery(j),UAVDelivery(k))~=inf
                combine=[combine;UAVDelivery(j),UAVDelivery(k)];
            end
            
        end
    end
    uncombine=setdiff(UAVDelivery,unique(combine)); % 不可合并的需求
    
    % 计算物资需求无人机的配送路径
    for j=1:size(combine,1)
        % 计算物资需求合并后的两个地点与车辆最优配送路径上的各个地点的距离
        d1=distance_min_2(path,combine(j,1));
        d2=distance_min_2(path,combine(j,2));
        
        % 计算物资需求合并后的两个地点与车辆途径地点的距离最近的地点，作为无人机起飞和降落地点
        D1=d1;D2=d2;
        e=[];
        e(1)=find(D1==min(D1),1);
        e(2)=find(D2==min(D2),1);
        % 无人机起飞和降落地点的两个地点重合，则进行调整
        if e(1)==e(2) && min(D1)<min(D2)
            D2(D2==min(D2))=max(D2);
            e(2)=find(D2==min(D2),1);
        elseif e(1)==e(2) && min(D1)>=min(D2)
            D1(D1==min(D1))=max(D1);
            e(1)=find(D1==min(D1),1);
        end
        path_UAV_mid(find(combine(j,1)==UAVDelivery,1),2)=path(min(e));
        path_UAV_mid(find(combine(j,1)==UAVDelivery,1),3)=path(max(e));
        path_UAV_mid(find(combine(j,2)==UAVDelivery,1),2)=path(min(e));
        path_UAV_mid(find(combine(j,2)==UAVDelivery,1),3)=path(max(e));
        
        % 计算无人机具体怎么走
        if distance_min_2(path(min(e)),combine(j,1))+distance_min_2(path(max(e)),combine(j,2))>=...
                distance_min_2(path(max(e)),combine(j,1))+distance_min_2(path(min(e)),combine(j,2))
            path_mid(find(UAVDelivery==combine(j,1),1),2:5)=[path(min(e)),combine(j,2),combine(j,1),path(max(e))]; % 记录第j个无人机配送最佳路径
            path_mid(find(UAVDelivery==combine(j,2),1),2:5)=[path(min(e)),combine(j,2),combine(j,1),path(max(e))]; % 记录第j个无人机配送最佳路径
        else
            path_mid(find(UAVDelivery==combine(j,1),1),2:5)=[path(min(e)),combine(j,1),combine(j,2),path(max(e))]; % 记录第j个无人机配送最佳路径
            path_mid(find(UAVDelivery==combine(j,2),1),2:5)=[path(min(e)),combine(j,1),combine(j,2),path(max(e))]; % 记录第j个无人机配送最佳路径
            
        end
    end
    
    % 计算不可合并需求的无人机配送路径
    d=[];
    for p=1:length(path)-1
        d=[d;distance_min_2(path(p),uncombine)+distance_min_2(path(p+1),uncombine)];
    end
    
    % 判断是否还可以交换，得到更优解
    e=[];
    for q=1:length(uncombine)
        e(q)=find(d(:,q)==min(d(:,q)),1);
        if q>=2 && ~isempty(find(e(q)==e(1:q-1),1)) % 有选位冲突
            D=d;
            E=e;
            SORT=sort(D(:,find(e(q)==e,1)));
            max_1=SORT(1);
            big_1=SORT(2);
            SORT=sort(D(:,q));
            max_2=SORT(1);
            big_2=SORT(2);
            if max_1+big_2>big_1+max_2 % 选位互换总和更小
                E(find(e(q)==e,1))=find(d(:,find(e(q)==e,1))==big_1,1);
                
                path_UAV_mid(find(uncombine(find(e(q)==e,1))==UAVDelivery,1),2)=path(E(find(E(q)==e,1)));
                path_UAV_mid(find(uncombine(find(e(q)==e,1))==UAVDelivery,1),3)=path(E(find(E(q)==e,1))+1);
                path_mid(find(uncombine(find(e(q)==e,1))==UAVDelivery,1),2:4)=path_UAV_mid(find(uncombine(find(e(q)==e,1))==UAVDelivery,1),[2,1,3]); % 记录无人机配送最佳路径
                path_UAV_mid(find(uncombine(q)==UAVDelivery,1),2)=path(E(q));
                path_UAV_mid(find(uncombine(q)==UAVDelivery,1),3)=path(E(q)+1);
                path_mid(find(uncombine(q)==UAVDelivery,1),2:4)=path_UAV_mid(find(uncombine(q)==UAVDelivery,1),[2,1,3]); % 记录无人机配送最佳路径
                e=E;
            else % 重新选位
                D=d;
                D(find(D(:,q)==min(D(:,q)),1),q)=max(D(:,q));
                e(q)=find(D(:,q)==min(D(:,q)),1);
                path_UAV_mid(find(uncombine(q)==UAVDelivery,1),2)=path(e(q));
                path_UAV_mid(find(uncombine(q)==UAVDelivery,1),3)=path(e(q)+1);
                path_mid(find(uncombine(q)==UAVDelivery,1),2:4)=path_UAV_mid(find(uncombine(q)==UAVDelivery,1),[2,1,3]); % 记录无人机配送最佳路径
            end
            
        else
            path_UAV_mid(find(uncombine(q)==UAVDelivery,1),2)=path(e(q));
            path_UAV_mid(find(uncombine(q)==UAVDelivery,1),3)=path(e(q)+1);
            path_mid(find(UAVDelivery==uncombine(q),1),2:4)=path_UAV_mid(find(UAVDelivery==uncombine(q),1),[2,1,3]); % 记录无人机配送最佳路径
            
        end
        
    end
    
    % 记录无人机最优配送方案
    path_UAV{i}=path_mid;
    
    
    
    
    % 如能完成，车辆额外等待的时间=无人机在两地点之间飞行时间-车辆最佳配送路线上的相同两地点之间行驶所花费的时间
    % 例如：汽车最佳配送路线中包含6-10，里程42km，耗时0.84h；无人机配送路线包含6-3-4-10，里程59，耗时0.7867h，此时车辆无需等待。
    VehicleWaitTime=0;
    % 如果属于combine任务，则只计算一次路径时间
    for j=1:size(combine,1)
        UAVDelivery(UAVDelivery==combine(j,2))=[];
    end
    % 对已经去掉合并需求的UAVDelivery进行计算
    for j=1:length(UAVDelivery)
        % 无人机在两地点之间飞行时间time1
        % 无人机配送任务的详细路径
        Path_UVA=path_mid(find(path_mid(:,1)==UAVDelivery(j),1),2:end);
        Path_UVA(Path_UVA==0)=[];
        % 无人机飞行路径距离计算
        Dist_UVA=0;
        for k=1:length(Path_UVA)-1
            Dist_UVA=Dist_UVA+distance_min_2(Path_UVA(k),Path_UVA(k+1));
        end
        % 无人机飞行时间
        time1=Dist_UVA/UAVSpeed;
        
        % 计算车辆从无人机起飞地点到降落地点的时间
        time2=0;
        Path_Vehicle=path(find(Path_UVA(1)==path, 1 ):find(Path_UVA(end)==path, 1, 'last' ));
        % 车辆配送路径距离计算
        Dist_Vehicle=0;
        for m=1:length(Path_Vehicle)-1
            Dist_Vehicle=Dist_Vehicle+distance_min_1(Path_Vehicle(m),Path_Vehicle(m+1));
        end
        % 车辆配送时间
        time2=Dist_Vehicle/VehicleSpeed;
        
        % 判断车辆是否需要等待，计算等待时间
        time=time1-time2;
        if time>0
            VehicleWaitTime=VehicleWaitTime+time;
        end
    end
    VehicleWaitTime;
    
    % 总时间=车辆最佳配送路线所耗时间+车辆额外等待的时间
    Time_Total=time_total+VehicleWaitTime;
    
    % 最短总时间即为目标函数
    objvalue(i)=Time_Total;
    
    
    
end

end