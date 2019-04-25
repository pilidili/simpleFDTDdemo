 %function Two_PML=Two_PML(time,freq)
time=3000;
freq=500;
t1=clock;

delta_x=0.01;
delta_y=0.01;
delta_t = delta_x/(2*3e8);
max_time = time;
max_space_x=400;
max_space_y=200;
max_space_x_si=200;
max_space_y_si=200;
%set up parameters for Si and SiO2
pi=3.1415926;
eps0=8.85419e-12;
mu0=4*pi*1e-7;          %磁导率
sigma_si=2.673e-4;      %电导率
eps_si=11.9*eps0;       %SiO2介电常数
eps_sio2=3.9*eps0;
%%初始化电磁场
Ez=zeros(max_space_x,max_space_y);
Hx=Ez;
Hy=Ez;
M=Ez;
%初始化常数
sigma=zeros(max_space_x,max_space_y);           %初始化电导率矩阵
for i=1:max_space_x
    for j=1:max_space_y
        if i<max_space_x_si && j<=max_space_y_si
            sigma(i,j)=sigma_si;
        else
            if i==max_space_x_si
                sigma(i,j)=sigma_si/2;          %界面取均值
            end
        end
    end
end

eps=zeros(max_space_x,max_space_y);             %初始化磁导率矩阵
for i=1:max_space_x
    for j=1:max_space_y
        if i<max_space_x_si &&j<=max_space_y_si
            eps(i,j)=eps_si;
        else
            if i==max_space_x_si
                eps(i,j)=(eps_si+eps_sio2)/2;   %边界取均值
            else
                eps(i,j)=eps_sio2;
            end
        end
    end
end

%初始化差分系数CA,CB,CP,CQ
CA=zeros(max_space_x,max_space_y);
for i=1:max_space_x
    for j=1:max_space_y
        CA(i,j)=(1-sigma(i,j)*delta_t/(2*eps(i,j)))/(1+sigma(i,j)*delta_t/(2*eps(i,j)));
    end
end
CB=zeros(max_space_x,max_space_y);
for i=1:max_space_x
    for j=1:max_space_y
        CB(i,j)=delta_t/eps(i,j)/(1+sigma(i,j)*delta_t/2/eps(i,j));
    end
end
CP=zeros(max_space_x,max_space_y);
for i=1:max_space_x
    for j=1:max_space_y
        CP(i,j)=(1-sigma(i,j)*delta_t/2/mu0)/(1+sigma(i,j)*delta_t/2/mu0);
    end
end
CQ=zeros(max_space_x,max_space_y);
for i=1:max_space_x
    for j=1:max_space_y
        CQ(i,j)=delta_t/mu0/(1+sigma(i,j)*delta_t/2/mu0);
    end
end
%初始化激励源
frequency=freq*1e6;                         % 频率单位为MHz
pulse_x=50;                                 %激励源位置
pulse_y=50;


figure(1)
set(gcf,'Position',get(0,'ScreenSize'));
ColorMap=load('colormap.txt');
%FDTD 循环开始
for T=1:max_time
    %加激励源
    pulse = sin(2*pi*frequency*delta_t*T);
    Ez(pulse_x,pulse_y)=Ez(pulse_x,pulse_y)+pulse;
    %更新H
    for i=1:max_space_x-1
        for j=1:max_space_y-1
            Hx(i,j)=CP(i,j)*Hx(i,j)-CQ(i,j)*(Ez(i,j+1)-Ez(i,j))/delta_y;
            Hy(i,j)=CP(i,j)*Hy(i,j)+CQ(i,j)*(Ez(i+1,j)-Ez(i,j))/delta_x;
        end
    end
    %更新E
    M=Ez;                                          %更新前保存上一次的Ez数据
    for i=2:(max_space_x-1)
        for j=2:(max_space_y-1)
            Ez(i,j)=CA(i,j)*Ez(i,j)+CB(i,j)*((Hy(i,j)-Hy(i-1,j))/delta_x-(Hx(i,j)-Hx(i,j-1))/delta_y);
        end
    end
    %Mur吸收边界（一阶近似）
    for j=2:(max_space_y-1)
        Ez(1,j)=M(2,j)-(Ez(2,j)-Ez(1,j))/3;        %左边界
        Ez(max_space_x,j)=M(max_space_x-1,j)-(Ez(max_space_x-1,j)-Ez(max_space_x,j))/3;    %右边界
    end
    for i=2:(max_space_x-1)
        Ez(i,1)=M(i,2)-(Ez(i,2)-Ez(i,1))/3;        %下边界
        Ez(i,max_space_y)=M(i,max_space_y-1)-(Ez(i,max_space_y-1)-Ez(i,max_space_y))/3;    %上边界
    end
    %二维角点处理
    Ez(1,1)=M(2,2)-0.47759*(Ez(2,2)-Ez(1,1));      %左下
    Ez(1,max_space_y)=M(2,max_space_y-1)-0.47759*(Ez(2,max_space_y-1)-Ez(1,max_space_y));%左上
    Ez(max_space_x,1)=M(max_space_x-1,2)-0.47759*(Ez(max_space_x-1,2)-Ez(max_space_x,1));%右下
    Ez(max_space_x,max_space_y)=M(max_space_x-1,max_space_y-1)-0.47759*(Ez(max_space_x-1,max_space_y-1)-Ez(max_space_x,max_space_y));
                                                                                         %右上
    %绘制Ez图像
    Ez_plot=Ez*100;
%     subplot(2,2,1)
%     surf(Ez_plot,'EdgeColor','None')
%     zlim([-100,100]);
%     view(0,0)
%     axis equal
%     axis([0,max_space_y,0,max_space_x,-100,100]);
    colormap(ColorMap);
    subplot(2,2,2)
    surf(Ez_plot,'EdgeColor','None')
    zlim([-100,100]);
    caxis([-80,80]);
    view(-90,0)
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
    
    subplot(2,2,4)
    surf(Ez_plot,'EdgeColor','None')
    zlim([-100,100]);
    caxis([-80,80]);
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
    view(-90,90)
    
    subplot(2,2,[1,3])
    surf(Ez_plot,'EdgeColor','None')
    zlim([-100,100]);
    caxis([-80,80]);
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
    t2=clock;
	totletime=etime(t2,t1);
    title(['程序运行：',num2str(totletime),' s，迭代次数:',num2str(T),'，平均用时:',num2str(totletime/T,2),' s/次']);

    pause(1e-20)
end











