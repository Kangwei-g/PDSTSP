% 2.3 ����������Ӧֵ
%�Ŵ��㷨�ӳ���
%Name:calfitvalue.m
%����������Ӧֵ
function fitvalue=ycsf_calfitvalue(objvalue,Tmean)

[px,py]=size(objvalue);
for i=1:px
    if Tmean-objvalue(i)>=0
        temp=Tmean-objvalue(i)+0.01;
    else
        temp=0.0;
    end
    fitvalue(i)=temp;
end
fitvalue=fitvalue';
