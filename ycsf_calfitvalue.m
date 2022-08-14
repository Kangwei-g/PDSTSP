% 2.3 计算个体的适应值
%遗传算法子程序
%Name:calfitvalue.m
%计算个体的适应值
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
