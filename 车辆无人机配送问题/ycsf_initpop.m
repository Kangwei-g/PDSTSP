% 2.1��ʼ��(����)
% initpop.m�����Ĺ�����ʵ��Ⱥ��ĳ�ʼ����popsize��ʾȺ��Ĵ�С
%�Ŵ��㷨�ӳ���
%Name: initpop.m
%��ʼ��
function pop=ycsf_initpop(popsize,chromlength)
% �Ŵ��㷨���룬�Զ��������ֳ��ȱ�ʾ���˻��������͵ĵص㣬0��ʾ�����͸õص㣬1��ʾ����
pop=[];
for t=1:popsize
    ChromosomeCode=ones(1,chromlength);
    site=round(rand(1)*chromlength);
    if site>0
        ChromosomeCode(site)=round((rand(1)+0.6)/2);
    end
    pop=[pop;ChromosomeCode];
    % rand�������ÿ����ԪΪ{0,1}
    % ����Ϊpopsize������Ϊchromlength�ľ���
    % round�Ծ����ÿ����Ԫ����Բ�������������ĳ�ʼ��Ⱥ��
end

end