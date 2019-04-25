%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright: Copyright (c) 2018
%Created on 2018-4-5  
%Author:dilipili
%Version 1.0 
%Title: FDTD_PML_MatrixVersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

time=3000;
freq=500;
t1=clock;%��ʼ��ʱ

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
mu0=4*pi*1e-7;          %�ŵ���
sigma_si=2.673e-4;      %�絼��
eps_si=11.9*eps0;       %SiO2��糣��
eps_sio2=3.9*eps0;
%%��ʼ����ų�
Ez=zeros(max_space_x,max_space_y);
Hx=Ez;
Hy=Ez;
M=Ez;
%��ʼ������
sigma=zeros(max_space_x,max_space_y);           %��ʼ���絼�ʾ���
sigma(1:max_space_x_si-1,1:max_space_y_si)=sigma_si;
sigma(max_space_x_si,:)=sigma_si/2;

eps=zeros(max_space_x,max_space_y);             %��ʼ���ŵ��ʾ���
eps(1:max_space_x_si-1,1:max_space_y_si)=eps_si;
eps(max_space_x_si,:)=(eps_si+eps_sio2)/2;
eps(max_space_x_si+1:max_space_x,:)=eps_sio2;

%��ʼ�����ϵ��CA,CB,CP,CQ
CA=(1-sigma.*delta_t./(2*eps))./(1+sigma.*delta_t./(2*eps));
CB=delta_t./eps./(1+sigma.*delta_t./2./eps);
CP=(1-sigma.*delta_t/2/mu0)./(1+sigma.*delta_t/2/mu0);
CQ=delta_t/mu0./(1+sigma.*delta_t/2/mu0);

%��ʼ������Դ
frequency=freq*1e6;                         % Ƶ�ʵ�λΪMHz
pulse_x=50;                                 %����Դλ��
pulse_y=50;


set(gcf,'Position',get(0,'ScreenSize'));
colormap(load('colormap.txt'));
%FDTD ѭ����ʼ
for T=1:max_time
    %�Ӽ���Դ
    pulse = sin(2*pi*frequency*delta_t*T);
    Ez(pulse_x,pulse_y)=Ez(pulse_x,pulse_y)+pulse;
    %����H
    i=1:max_space_x-1;
    j=1:max_space_y-1;
    Hx(i,j)=CP(i,j).*Hx(i,j)-CQ(i,j).*(Ez(i,j+1)-Ez(i,j))/delta_y;
    Hy(i,j)=CP(i,j).*Hy(i,j)+CQ(i,j).*(Ez(i+1,j)-Ez(i,j))/delta_x;
    clear i j;
    
    %����E
    M=Ez;                                          %����ǰ������һ�ε�Ez����
    i=2:(max_space_x-1);
    j=2:(max_space_y-1);
    Ez(i,j)=CA(i,j).*Ez(i,j)+CB(i,j).*((Hy(i,j)-Hy(i-1,j))/delta_x-(Hx(i,j)-Hx(i,j-1))/delta_y);
    clear i j;

    %Mur���ձ߽磨һ�׽��ƣ�
    j=2:(max_space_y-1);
    Ez(1,j)=M(2,j)-(Ez(2,j)-Ez(1,j))/3;        %��߽�
    Ez(max_space_x,j)=M(max_space_x-1,j)-(Ez(max_space_x-1,j)-Ez(max_space_x,j))/3;    %�ұ߽�
    
    i=2:(max_space_x-1);
    Ez(i,1)=M(i,2)-(Ez(i,2)-Ez(i,1))/3;        %�±߽�
    Ez(i,max_space_y)=M(i,max_space_y-1)-(Ez(i,max_space_y-1)-Ez(i,max_space_y))/3;    %�ϱ߽�


%��ά�ǵ㴦��
    Ez(1,1)=M(2,2)-0.47759*(Ez(2,2)-Ez(1,1));      %����
    Ez(1,max_space_y)=M(2,max_space_y-1)-0.47759*(Ez(2,max_space_y-1)-Ez(1,max_space_y));%����
    Ez(max_space_x,1)=M(max_space_x-1,2)-0.47759*(Ez(max_space_x-1,2)-Ez(max_space_x,1));%����
    Ez(max_space_x,max_space_y)=M(max_space_x-1,max_space_y-1)-0.47759*(Ez(max_space_x-1,max_space_y-1)-Ez(max_space_x,max_space_y));
                                                                                         %����
    %����Ezͼ��
    Ez_plot=Ez*100;
    subplot(2,2,2)
    surf(Ez_plot,'EdgeColor','None')
    line([1,1],[max_space_x_si,max_space_x_si],[-100,100]);
    zlim([-100,100]);
    caxis([-80,80]);
    view(-90,0)
    title('����ͼ')
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
        
    subplot(2,2,4)
    surf(Ez_plot,'EdgeColor','None')
    line([1,max_space_y],[max_space_x_si,max_space_x_si]);
    zlim([-100,100]);
    caxis([-80,80]);
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
    view(-90,90)
    title('����ͼ')
    
    subplot(2,2,[1,3])
    surf(Ez_plot,'EdgeColor','None')
    line([1,max_space_y],[max_space_x_si,max_space_x_si]);
    zlim([-100,100]);
    caxis([-80,80]);
    axis equal
    axis([0,max_space_y,0,max_space_x,-100,100]);
    t2=clock;
	totletime=etime(t2,t1);
    title(['�������У�',num2str(totletime),' s������������',num2str(T),'��ƽ��',num2str(totletime/T,2),' s/�Σ�ʣ��ʱ�䣺',num2str(totletime/T*(time-T),3),' s']);

    pause(1e-20)
end
    subplot(2,2,[1,3])
    title(['�������У�',num2str(totletime),' s������������',num2str(T),'��ƽ��',num2str(totletime/T,2),' s/�Σ�������ɣ�']);







