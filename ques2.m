clc,clear
%% 2022��繤��B�������
% ͼ 2 ��ʵ�ߴ����������˻��������ߵ�·�ߣ����ߴ���ֻ�����˻������ߵ�·�ߡ�
% Ӧ��������Ȼ�����ڵ� 9 ���ص㣬���ͳ��������������Ϊ 1000 ǧ�ˣ���ȡ�����ͳ���+���˻���������ģʽ��
%% Ŀ�꺯�������һ��������������Ҫ��ʱ�����
%% ���߱���
% ���ʼ��еص�
% ����Ŀ��ص�
% ����֮��·�����
% ����֮�����
% ���㵱������������
% ���ͳ������������

%% Լ������
% ���ͳ���������
% ���ͳ�������ʻ�ٶ�
% ���˻�������ÿ�����ʼ��е�һ̨��
% ���˻����������
% ���˻���ʻ�ٶ�
% ���˻������ʱ��

%% ��һ���������ݽ��д�������matlab
%ÿ������������
MaterialNeed=[12 90 24 15 70 18 150 50 30 168 36 44 42 13]; % ǰ���ʵ�����

% ��Ҫ���͵ĵص�
Position=[1:1:14];

% ��Ҫ���͵ĵص�����
PositionNum=length(Position);

% Ӧ�����ʼ��еص�
DistributionCenter=9;

% ���ʼ��еص�����
DistributionCenterNum=length(DistributionCenter);

% �����֮�����
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
   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    36    25   Inf   Inf    47   Inf];  % �����ľ��루���������˻�·����

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
   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf    36    25   Inf   Inf    47   Inf];  % ���ӵ����˻���·����Ҳ�����˻����Է��е�·��


% ���ͳ������������
VehicleCapacity=1000; % ����1��2��������������

% ���ͳ�������
VehicleNum=DistributionCenterNum; 

% ���ͳ�����ʻ�ٶ�km/h

VehicleSpeed=50; 

% ���˻���������������˻�unmanned aerial vehicles��UAV

UAVCapacity=50; 
% ���˻�����
UAVNum=VehicleNum; 
% ���˻���ʻ�ٶ�km/h
UAVSpeed=75;
% ���˻������ʱ��h
UAVEnduranceTime=70/60;
% ���˻���󺽳�km
UAVEnduranceDistance=UAVSpeed*UAVEnduranceTime;


% �����Ѿ�����õ�����
load('data_of_2022527B_ques2.mat')

%% �ڶ������ȼ��������ص�֮�����̾��롢���·�������ʱ��
% ���ͳ�������ʻ·��
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

% ���˻�����ʻ·�����˴����ǵ����ʱ��ʵ����Ҫ�����ǻ���Ҫ���������˻���������ڵ����пɷ���·��
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

% ƴ������path_min_2Ϊһ��
record_path_min_2={};
record_distance_min_2=[];
for i=1:PositionNum
        record_path_min_2=[record_path_min_2;path_min_2(:,i)];
        record_distance_min_2=[record_distance_min_2;distance_min_2(:,i)];
end

% ���������˻���󺽳̵�·����Ϊ�����ֻ�������˻�����������·��
a=find(record_distance_min_2>UAVEnduranceDistance);
record_path_min_2(a,:)=[];
record_distance_min_2(a,:)=[];

% ���������˻���������ڵ����пɷ���·��
path_temp=record_path_min_2;
path_end={};
distance_temp=record_distance_min_2;
distance_end=[];
while ~isempty(path_temp) % ��ʱ��¼δ����
    % ��Զ�ӵ�һ��·�������ж�
    pos_last=path_temp{1}(end); % ��ǰ·�������ص�
    b=find(Distance_2(pos_last,:)~=inf);
    distance_sum=[];
    for i=1:length(b)
        distance_new=distance_temp(1)+Distance_2(pos_last,b(i));
        distance_sum=[distance_sum;distance_new];
        if distance_new<=UAVEnduranceDistance % ��·���ľ���<=���˻�����������Լ�����¼·��
            distance_temp=[distance_temp;distance_new];
            path_new=[path_temp{1},b(i)];
            path_temp=[path_temp;path_new];
            
            % ��·���ľ���<=���˻����������·��Ϊһ������·�������м�¼
            path_end=[path_end;path_new];
            distance_end=[distance_end;distance_new];
        end
    end

    % ��һ��·���Ѿ�������ϣ�ɾ��
    path_temp(1)=[];
    distance_temp(1)=[];
end

% ��path_end��Ԫ��ת��Ϊ����
path_end_array=[];
for i=1:length(path_end)
    path_end_array(i,1:length(path_end{i}))=path_end{i};
end
path_end=path_end_array;

%% ��������������Щ�ص���������˻�����
% ��Щ�ص���������˻����ͣ����Ǽ������˻�����������ͬһ���ص�����
% ��Ϊ�����ͬһ���ط����ͣ�������ʻ�ٶ�50km/h�����˻������ٶ�75km/h���������εĻ��������ӵĵȴ�ʱ��>�����˻����ͽ�ʡ��ʱ��
% ����ÿ������������>���˻����������50kg�ĵص㶼�������˻����͵Ŀ��Ƿ�Χ�ڡ�
% ���ܵ����˻����͵ص�

UAVDeliveryPositionMaybe=find(MaterialNeed<=UAVCapacity);
UAVDeliveryPositionMaybe(UAVDeliveryPositionMaybe==DistributionCenter)=[]; % ���ʼ��еز��������˻����ͣ�����ܵ����˻����͵ص������ʼ��еأ�����ɾ��

% ���˻��������͵ĵص�

UAVNoFlyPosition=find(MaterialNeed>UAVCapacity);
UAVNoFlyPosition=sort(unique([UAVNoFlyPosition,DistributionCenter])); % ���ʼ��еز��������˻����ͣ������˻��������͵ĵص�û�����ʼ��еأ��������


save ques2_calculated_data.mat MaterialNeed Distance_1 Distance_2 VehicleSpeed UAVCapacity UAVSpeed distance_min_1 path_min_1 time_min_1 ...
    distance_min_2 path_min_2 time_min_2 path_end distance_end UAVDeliveryPositionMaybe UAVNoFlyPosition
%% ���Ĳ��������˻����͵ĵص�ΪȾɫ�壬ʹ���Ŵ��㷨���ĳһ�����Ŀ�꺯���������һ��������������Ҫ��ʱ�����

popsize=100; %Ⱥ���С
chromlength=length(UAVDeliveryPositionMaybe); %�ַ������ȣ����峤��Ϊ���˻��������͵ĵص�����
pc=0.001; %�������
pm=0.0001; %�������
Tmean=10; % ��Ӧ�ȣ�����ʱ�䣩


chromlength=length(UAVDeliveryPositionMaybe); % ���˻������͵ص�����ΪȾɫ�峤��
pop=ycsf_initpop(popsize,chromlength); %������ʼȺ��

for i=1:5 %��������
    [objvalue,~,~,~]=ycsf_calobjvalue(pop); %����Ŀ�꺯��
    Tmean=mean(objvalue); % ��Ӧ�ȣ�����ȡƽ������ʱ��
    fitvalue=ycsf_calfitvalue(objvalue,Tmean); %����Ⱥ����ÿ���������Ӧ��
    [newpop]=ycsf_selection(pop,fitvalue);  % ѡ��
    [newpop]=ycsf_crossover(pop,pc);  % ����
    [newpop]=ycsf_mutation(pop,pm);   % ����
    pop=newpop;
end
[bestindividual,bestfit]=ycsf_best(pop,fitvalue); %���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ

%% ���岽���������Ѹ�������ų������ͷ��������˻����ͷ������������ͷ���������̡��������ʱ��
[objvalue,dist_vehicle,path_vehicle,path_UAV]=ycsf_calobjvalue(bestindividual);
while isempty(path_UAV)
    [objvalue,dist_vehicle,path_vehicle,path_UAV]=ycsf_calobjvalue(bestindividual);
end
%% չʾ�ͱ�������

disp('�����������ͷ���path_vehicle=')
path_vehicle{1}
disp('���˻��������ͷ�������һ�������������˻����͵ĵص㣬֮��������˻���ɡ�;���ͽ���ĵص㣩path_UAV=')
path_UAV{1:end}
disp('�����������ͷ����������dist_vehicle=')
dist_vehicle
disp('�������ʱ��MinTime=')
MinTime=objvalue

save ques2_result.mat path_vehicle path_UAV dist_vehicle MinTime 
