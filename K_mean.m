function [center,guilei]=K_mean(k,Distance_1,distance_min_1)

k ; %��������


point_num=size(Distance_1,1);
point_num; % dΪ������ĸ���

%��һ������ʼ���ѡ��k������Ϊ���ġ�
%��һ���������������һ���⣬��Ϊ��ʼ���н⡣�ʳ�ʼ���Ӱ�������ٶȺ����վ�����������ֲ����Ž⣩������ͨ��������ѡȡ��ʼ���ҵ����ŵĽ����

center = zeros(k,1);
suiji = randperm(point_num);

% center = suiji(1:k)';
center = [6,7];


ddcs = 50; %���õ�������������������һ��������ֹͣ
zjl = zeros(1,ddcs); %�ܾ��룬�������ĵ����ڸ������֮�͡�

for q = 1:1:ddcs
    %�ڶ����������е㣬������ֵ������������õ�k���ࡣ
    %�ڶ�����������������������������ġ�
    
    %����ÿ���㵽���ĵľ���
    distences = zeros(k,point_num);
    dist=distance_min_1;
    dist(dist==inf)=0;
    for i = 1:1:k
	for j = 1:1:point_num
		distences(i,j) = dist(center(i),j);
	end
    end
    %����������̾���ĵ����
    guilei = zeros(1,point_num);
    for i = 1:1:point_num
        guilei(1,i) = min(find(distences(:,i) == min(distences(:,i))));
    end
    
    
    
    
    %������������Ѱ��������
    %����������������һƬ���ݵ���˵�����ݵ�����ģ������ƽ��ֵ��������Щ����ܾ�����̡�

    %�����ĵ�ԭ���������е�������  
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
    %���Ĳ���������û�иı䣬��������
    %�������������ڵ���ܾ���
    %���Ĳ������������ĸı䣬���ܾ���仯�����ܾ���ı仯���ж������Ƿ��б仯��
    zjl(1,q) = sum(sum(DistSum));
    if q > 1
        if zjl(1,q) == zjl(1,q - 1)
            break
        end
    end
    
end