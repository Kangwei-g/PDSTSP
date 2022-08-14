clc,clear
%% 2022年电工杯B题问题二
% 图 2 中实线代表车辆和无人机都可以走的路线，虚线代表只有无人机可以走的路线。
% 应急物资仍然集中在第 9 个地点，配送车辆的最大载重量为 1000 千克，采取“配送车辆+无人机”的配送模式。
%% 目标函数：完成一次整体配送所需要的时间最短
%% 决策变量
% 物资集中地点
% 配送目标地点
% 各点之间路线情况
% 各点之间距离
% 各点当日物资需求量
% 配送车辆最大载重量

%% 约束条件
% 配送车辆的数量
% 配送车辆的行驶速度
% 无人机数量（每个物资集中点一台）
% 无人机最大载重量
% 无人机行驶速度
% 无人机最长飞行时间

%% 第一步：将数据进行处理，导入matlab
%每日物资需求量
MaterialNeed=[12 90 24 15 70 18 150 50 30 168 36 44 42 13]; % 前三问的数据

% 需要配送的地点
Position=[1:1:14];

% 需要配送的地点数量
PositionNum=length(Position);

% 应急物资集中地点
DistributionCenter=9;

% 物资集中地点数量
DistributionCenterNum=length(DistributionCenter);

% 点与点之间距离
Distance_1=[Inf   Inf   Inf   Inf    54   Inf    55   Inf   Inf   Inf    26   Inf   Inf   Inf
   Inf   Inf    56   Inf    18   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf
   Inf    56   Inf   Inf    44   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf
   Inf   Inf   Inf   Inf   Inf    28   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf
    54    18    44   Inf   Inf    51    34    56    48   Inf   Inf   Inf   Inf   Inf
   Inf   Inf   Inf    28    51   Inf   Inf   Inf    27    42   Inf   Inf   Inf   Inf
    55   Inf   Inf   Inf    34   Inf   Inf    36   Inf   Inf   Inf    38   Inf   Inf
   Inf   Inf   Inf   Inf    56   Inf    36   Inf    29   Inf   Inf    33   Inf   Inf
   Inf   Inf   Inf   Inf    48    27   Inf    29   Inf    61   Inf    29    42    36
   Inf   Inf   Inf   Inf   Inf    42   Inf   Inf    61   Inf   Inf   Inf   Inf    25
    26   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    24   Inf   Inf
   Inf   Inf   Inf   Inf   Inf   Inf    38    33    29   Inf    24   Inf   Inf   Inf
   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    42   Inf   Inf   Inf   Inf    47
   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    36    25   Inf   Inf    47   Inf];  % 点与点的距离（不包括无人机路径）

Distance_2=[Inf    20   Inf   Inf    54   Inf    55   Inf   Inf   Inf    26   Inf   Inf   Inf
    20   Inf    56   Inf    18   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf
   Inf    56   Inf    15    44    18   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf
   Inf   Inf    15   Inf   Inf    28   Inf   Inf   Inf    26   Inf   Inf   Inf   Inf
    54    18    44   Inf   Inf    51    34    56    48   Inf   Inf   Inf   Inf   Inf
   Inf   Inf    18    28    51   Inf   Inf   Inf    27    42   Inf   Inf   Inf   Inf
    55   Inf   Inf   Inf    34   Inf   Inf    36   Inf   Inf    26    38   Inf   Inf
   Inf   Inf   Inf   Inf    56   Inf    36   Inf    29   Inf   Inf    33    25   Inf
   Inf   Inf   Inf   Inf    48    27   Inf    29   Inf    61   Inf    29    42    36
   Inf   Inf   Inf    26   Inf    42   Inf   Inf    61   Inf   Inf   Inf   Inf    25
    26   Inf   Inf   Inf   Inf   Inf    26   Inf   Inf   Inf   Inf    24   Inf   Inf
   Inf   Inf   Inf   Inf   Inf   Inf    38    33    29   Inf    24   Inf    34   Inf
   Inf   Inf   Inf   Inf   Inf   Inf   Inf    25    42   Inf   Inf    34   Inf    47
   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    36    25   Inf   Inf    47   Inf];  % 增加的无人机的路径，也即无人机可以飞行的路径


% 配送车辆最大载重量
VehicleCapacity=1000; % 问题1、2中汽车载重能力

% 配送车辆数量
VehicleNum=DistributionCenterNum; 

% 配送车辆行驶速度km/h

VehicleSpeed=50; 

% 无人机最大载重量，无人机unmanned aerial vehicles，UAV

UAVCapacity=50; 
% 无人机数量
UAVNum=VehicleNum; 
% 无人机行驶速度km/h
UAVSpeed=75;
% 无人机最长飞行时间h
UAVEnduranceTime=70/60;
% 无人机最大航程km
UAVEnduranceDistance=UAVSpeed*UAVEnduranceTime;


% 加载已经处理好的数据
load('data_of_2022527B_ques2.mat')

%% 第二步：先计算两两地点之间的最短距离、最短路径和最短时间
% 配送车辆可行驶路径
distance_min_1=[];
path_min_1={};
time_min_1=[];
for i=1:PositionNum
    for j=1:PositionNum
        [dist,path]=myfloyd(Distance_1,i,j);
        distance_min_1(i,j)=dist;
        time_min_1(i,j)=dist/VehicleSpeed;
        if j>=i
        path_min_1{i,j}=path;
        path_min_1{j,i}=flip(path);
        end
    end
end
distance_min_1;
path_min_1;
time_min_1;

% 无人机可行驶路径（此处考虑到求解时的实际需要，我们还需要计算在无人机续航里程内的所有可飞行路径
distance_min_2=[];
path_min_2={};
for i=1:PositionNum
    for j=1:PositionNum
        [dist,path]=myfloyd(Distance_2,i,j);
        distance_min_2(i,j)=dist;
        time_min_2(i,j)=dist/UAVSpeed;
        if j>=i
            path_min_2{i,j}=path;
            path_min_2{j,i}=flip(path);
        end
    end
end
distance_min_2;
path_min_2;
time_min_2;

% 拼接所有path_min_2为一列
record_path_min_2={};
record_distance_min_2=[];
for i=1:PositionNum
        record_path_min_2=[record_path_min_2;path_min_2(:,i)];
        record_distance_min_2=[record_distance_min_2;distance_min_2(:,i)];
end

% 将大于无人机最大航程的路径设为无穷大，只留下无人机可以往返的路径
a=find(record_distance_min_2>UAVEnduranceDistance);
record_path_min_2(a,:)=[];
record_distance_min_2(a,:)=[];

% 计算在无人机续航里程内的所有可飞行路径
path_temp=record_path_min_2;
path_end={};
distance_temp=record_distance_min_2;
distance_end=[];
while ~isempty(path_temp) % 临时记录未清零
    % 永远从第一个路径进行判断
    pos_last=path_temp{1}(end); % 当前路径的最后地点
    b=find(Distance_2(pos_last,:)~=inf);
    distance_sum=[];
    for i=1:length(b)
        distance_new=distance_temp(1)+Distance_2(pos_last,b(i));
        distance_sum=[distance_sum;distance_new];
        if distance_new<=UAVEnduranceDistance % 新路径的距离<=无人机续航，则可以继续记录路径
            distance_temp=[distance_temp;distance_new];
            path_new=[path_temp{1},b(i)];
            path_temp=[path_temp;path_new];
            
            % 新路径的距离<=无人机续航，则此路径为一个最终路径，进行记录
            path_end=[path_end;path_new];
            distance_end=[distance_end;distance_new];
        end
    end

    % 第一个路径已经处理完毕，删掉
    path_temp(1)=[];
    distance_temp(1)=[];
end

% 将path_end由元胞转变为数组
path_end_array=[];
for i=1:length(path_end)
    path_end_array(i,1:length(path_end{i}))=path_end{i};
end
path_end=path_end_array;

%% 第三步：计算哪些地点可以由无人机配送
% 哪些地点可以由无人机配送，我们假设无人机不会连续向同一个地点配送
% 因为如果向同一个地方配送，车辆行驶速度50km/h，无人机飞行速度75km/h，配送两次的话车辆增加的等待时间>有无人机配送节省的时间
% 所以每日物资需求量>无人机最大载重量50kg的地点都不在无人机配送的考虑范围内。
% 可能的无人机配送地点

UAVDeliveryPositionMaybe=find(MaterialNeed<=UAVCapacity);
UAVDeliveryPositionMaybe(UAVDeliveryPositionMaybe==DistributionCenter)=[]; % 物资集中地不可由无人机配送，如可能的无人机配送地点有物资集中地，则将其删除

% 无人机不能配送的地点

UAVNoFlyPosition=find(MaterialNeed>UAVCapacity);
UAVNoFlyPosition=sort(unique([UAVNoFlyPosition,DistributionCenter])); % 物资集中地不可由无人机配送，如无人机不能配送的地点没有物资集中地，则将其添加


save ques2_calculated_data.mat MaterialNeed Distance_1 Distance_2 VehicleSpeed UAVCapacity UAVSpeed distance_min_1 path_min_1 time_min_1 ...
    distance_min_2 path_min_2 time_min_2 path_end distance_end UAVDeliveryPositionMaybe UAVNoFlyPosition
%% 第四步：以无人机配送的地点为染色体，使用遗传算法求解某一组合下目标函数，即完成一次整体配送所需要的时间最短

popsize=100; %群体大小
chromlength=length(UAVDeliveryPositionMaybe); %字符串长度（个体长度为无人机可以配送的地点数）
pc=0.001; %交叉概率
pm=0.0001; %变异概率
Tmean=10; % 适应度（配送时间）


chromlength=length(UAVDeliveryPositionMaybe); % 无人机可配送地点数量为染色体长度
pop=ycsf_initpop(popsize,chromlength); %产生初始群体

for i=1:5 %迭代次数
    [objvalue,~,~,~]=ycsf_calobjvalue(pop); %计算目标函数
    Tmean=mean(objvalue); % 适应度，这里取平均配送时间
    fitvalue=ycsf_calfitvalue(objvalue,Tmean); %计算群体中每个个体的适应度
    [newpop]=ycsf_selection(pop,fitvalue);  % 选择
    [newpop]=ycsf_crossover(pop,pc);  % 交叉
    [newpop]=ycsf_mutation(pop,pm);   % 变异
    pop=newpop;
end
[bestindividual,bestfit]=ycsf_best(pop,fitvalue); %求出群体中适应值最大的个体及其适应值

%% 第五步：计算出最佳个体的最优车辆配送方案、无人机配送方案、最优配送方案的总里程、最短配送时间
[objvalue,dist_vehicle,path_vehicle,path_UAV]=ycsf_calobjvalue(bestindividual);
while isempty(path_UAV)
    [objvalue,dist_vehicle,path_vehicle,path_UAV]=ycsf_calobjvalue(bestindividual);
end
%% 展示和保存数据

disp('车辆最优配送方案path_vehicle=')
path_vehicle{1}
disp('无人机最优配送方案（第一个数字是由无人机配送的地点，之后的是无人机起飞、途径和降落的地点）path_UAV=')
path_UAV{1:end}
disp('车辆最优配送方案的总里程dist_vehicle=')
dist_vehicle
disp('最短配送时间MinTime=')
MinTime=objvalue

save ques2_result.mat path_vehicle path_UAV dist_vehicle MinTime 
