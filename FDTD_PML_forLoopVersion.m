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
for i=1:max_space_x
    for j=1:max_space_y
        if i<max_space_x_si && j<=max_space_y_si
            sigma(i,j)=sigma_si;
        else
            if i==max_space_x_si
                sigma(i,j)=sigma_si/2;          %����ȡ��ֵ
            end
        end
    end
end

eps=zeros(max_space_x,max_space_y);             %��ʼ���ŵ��ʾ���
for i=1:max_space_x
    for j=1:max_space_y
        if i<max_space_x_si &&j<=max_space_y_si
            eps(i,j)=eps_si;
        else
            if i==max_space_x_si
                eps(i,j)=(eps_si+eps_sio2)/2;   %�߽�ȡ��ֵ
            else
                eps(i,j)=eps_sio2;
            end
        end
    end
end

%��ʼ�����ϵ��CA,CB,CP,CQ
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
%��ʼ������Դ
frequency=freq*1e6;                         % Ƶ�ʵ�λΪMHz
pulse_x=50;                                 %����Դλ��
pulse_y=50;


figure(1)
set(gcf,'Position',get(0,'ScreenSize'));
ColorMap=load('colormap.txt');
%FDTD ѭ����ʼ
for T=1:max_time
    %�Ӽ���Դ
    pulse = sin(2*pi*frequency*delta_t*T);
    Ez(pulse_x,pulse_y)=Ez(pulse_x,pulse_y)+pulse;
    %����H
    for i=1:max_space_x-1
        for j=1:max_space_y-1
            Hx(i,j)=CP(i,j)*Hx(i,j)-CQ(i,j)*(Ez(i,j+1)-Ez(i,j))/delta_y;
            Hy(i,j)=CP(i,j)*Hy(i,j)+CQ(i,j)*(Ez(i+1,j)-Ez(i,j))/delta_x;
        end
    end
    %����E
    M=Ez;                                          %����ǰ������һ�ε�Ez����
    for i=2:(max_space_x-1)
        for j=2:(max_space_y-1)
            Ez(i,j)=CA(i,j)*Ez(i,j)+CB(i,j)*((Hy(i,j)-Hy(i-1,j))/delta_x-(Hx(i,j)-Hx(i,j-1))/delta_y);
        end
    end
    %Mur���ձ߽磨һ�׽��ƣ�
    for j=2:(max_space_y-1)
        Ez(1,j)=M(2,j)-(Ez(2,j)-Ez(1,j))/3;        %��߽�
        Ez(max_space_x,j)=M(max_space_x-1,j)-(Ez(max_space_x-1,j)-Ez(max_space_x,j))/3;    %�ұ߽�
    end
    for i=2:(max_space_x-1)
        Ez(i,1)=M(i,2)-(Ez(i,2)-Ez(i,1))/3;        %�±߽�
        Ez(i,max_space_y)=M(i,max_space_y-1)-(Ez(i,max_space_y-1)-Ez(i,max_space_y))/3;    %�ϱ߽�
    end
    %��ά�ǵ㴦��
    Ez(1,1)=M(2,2)-0.47759*(Ez(2,2)-Ez(1,1));      %����
    Ez(1,max_space_y)=M(2,max_space_y-1)-0.47759*(Ez(2,max_space_y-1)-Ez(1,max_space_y));%����
    Ez(max_space_x,1)=M(max_space_x-1,2)-0.47759*(Ez(max_space_x-1,2)-Ez(max_space_x,1));%����
    Ez(max_space_x,max_space_y)=M(max_space_x-1,max_space_y-1)-0.47759*(Ez(max_space_x-1,max_space_y-1)-Ez(max_space_x,max_space_y));
                                                                                         %����
    %����Ezͼ��
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
    title(['�������У�',num2str(totletime),' s����������:',num2str(T),'��ƽ����ʱ:',num2str(totletime/T,2),' s/��']);

    pause(1e-20)
end











