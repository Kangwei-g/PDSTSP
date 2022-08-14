% 2.1初始化(编码)
% initpop.m函数的功能是实现群体的初始化，popsize表示群体的大小
%遗传算法子程序
%Name: initpop.m
%初始化
function pop=ycsf_initpop(popsize,chromlength)
% 遗传算法编码，以二进制数字长度表示无人机可以配送的地点，0表示不配送该地点，1表示配送
pop=[];
for t=1:popsize
    ChromosomeCode=ones(1,chromlength);
    site=round(rand(1)*chromlength);
    if site>0
        ChromosomeCode(site)=round((rand(1)+0.6)/2);
    end
    pop=[pop;ChromosomeCode];
    % rand随机产生每个单元为{0,1}
    % 行数为popsize，列数为chromlength的矩阵，
    % round对矩阵的每个单元进行圆整。这样产生的初始种群。
end

end