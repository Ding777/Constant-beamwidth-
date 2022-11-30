clear
clc
close all;
warning off;
clear all; close all; clc;
%% 基本参数
f01 = 4000;           % 信号中心频率Signal center frequency

Nfft =512;                                % fft点数
F0=linspace(0, 8000, Nfft/2+1);           % The actual frequency corresponding to each fft frequency point
f_range=linspace(f01-B/2, f01+B/2, Nfft/2+1);     % The actual frequency corresponding to each fft frequency point
% F(1:10)=F(11);
% F(1:3)=F(4);

% F(1:16)=F(17);
% F(129:257)=F(128);
% F(1:end)=F(70);

% F(200:257)=F(199);
% F(1:128)=F(129);
% F = fl:B/(Nfft/2):fh; 

M =16;                                % Number of array elements
c = 340;                              % speed of sound
lambda = c/f01;                        %center wavelength
R =0.025;                               % actual radius
% R =0.5*M*lambda/(2*pi);             % 理想半径

gamma=2*pi*[0:M-1]'/M ;
% d = 0.5*lambda;                        % 阵元间距(半波长为理想宽度，可视距离内无额外峰值)
% BW= 0.886*lambda/(M*d)*180/pi;         % 中心频率处半功率波束宽度（分辨率）
% fprintf('半功率波束宽度BW = %0.2f°\n',BW)

% seta1=[60]*pi/180;                     % 零点约束俯仰角1
% fai1=[-140]*pi/180;                    % 零点约束方位角1  
 
% seta2=[60]*pi/180;                     % 零点约束俯仰角2
% fai2=[-111]*pi/180;                    % 零点约束方位角2 

%% 读入实际音频数据
seta0 = 60*pi/180;                       % Signal Arrival Elevation Angle
fai0  =-150*pi/180;                       % Signal Arrival Azimuth
f=[1]';    
%% diffuse noise field MSC
Fs=16000;                                                               
f0 = (1:Nfft/2+1)*Fs/Nfft;
Fvv = zeros(Nfft/2+1,M,M);
 %f0=500;
% d=2*R*sin(pi/M);

for i = 1:M
    for j = 1:M   
        if i == j
            Fvv(:,i,j) = ones(1,Nfft/2+1);
        else
            dij = 2*R*sin(pi*abs(i-j)/M);
            Fvv(:,i,j) = sin(2*pi*f0*dij*1.0/c)./(2*pi*f0*dij*1.0/c);
        end
    end
end

% 对角加载映射曲线
u=1:1:Nfft/2+1;
k=-0.3;
B=80;
Amin=0.1;
P=(1-Amin)./(1+exp(k*(u-B)))+Amin;     
% plot(u,P);ylim([-0.1,1])
for n=1:Nfft/2+1
%         Fvv2(n,:,:)=inv(squeeze(Fvv(n,:,:))+eye(M)*P(n));
      Fvv2(n,:,:)=inv(squeeze(Fvv(n,:,:))+eye(M)*0.1);
end
%  imshow(squeeze(Fvv(200,:,:)))
%  for i=1:256
%     imshow(squeeze(Fvv(i,:,:)))
%     pause(0.05)
%  end
%M = 10;                 %阵元数
%dd= (0:1:M-1) * (c/f0/2);   %阵元位置
p_all=[];
a0_all=[];
for nf =2:Nfft/2+1
    freq = f_range(nf);
    freq1=500;
   % a0 = exp(-1i * 2 * pi * freq /c *dd.'*sind(theta0));   %计算期望导向矢量
    a0=exp(-1j*2*pi*R*freq/c*sin(seta0)*cos(fai0-gamma));   %Calculate desired steering vector
    a0=a0;
    C=[a0];
   % a0=a0/100;
    Rxx=squeeze(Fvv2(nf,:,:));
%   Rxx=1;
    Wc=(Rxx*C)/(C'*Rxx*C)*f;                % MVDR algorithm to obtain weight coefficient
%       零点调向(支持多个零点)
    Wopt=Wc;  
    step = 1;
    seta=seta0; 
    scan_angle = [-180:step:180];                        %Scan Angle Range
    scan_angle=[-180:step:180]/180*pi; 
    a_all =[];                           
    a_all1 =[];  
    a_all2=[];
    for i_1 = 1:length(scan_angle) 
        %as = exp(-1i * 2 * pi * freq /c *dd.'*sind(scan_angle(i_1)));
        as=exp(-1j*2*pi*R*freq/c*sin(seta)*cos(scan_angle(i_1)-gamma')); 
        p1(i_1,nf)=abs(as*conj(Wopt));
       % p11(i_1,nf)=20*log10(abs(as*conj(Wopt))/max(abs(as*conj(Wopt))));
        %p11(nf,i_1)=abs(as);
       % p4(i_1,nf)=20*log10(abs(as*conj(Wopt))/max(abs(as*conj(Wopt))));
        a_all= [a_all as'];
        a_all1= [a_all1  (20*log10(abs(as)/max(abs(as))))]; 
        a_all2=[a_all2  (20*log10(abs(as.*conj(Wopt))/max(abs(as.*conj(Wopt)))))];
    end
    %p11(i_1,nf)=abs(a_all);
    p = 20*log10(abs(a_all'*conj(Wopt))/max(abs(a_all'*conj(Wopt))));  %Computational pattern
    p3 = 20*log10(abs(a_all')/max(abs(a_all')));  
    p4 = 20*log10(abs(p1)/max(abs(p1))); 
    %p = 20*log10(abs(a_all'*conj(a0))/max(abs(a_all'*conj(a0))));  
    p7 = 20*log10(abs(a_all'*(Wopt))/max(abs(a_all'*(Wopt))));
    p11=conj(a0)';
        f01=5.312500000000000e+03;
        f01=1000;
        if freq == f01
        angle_ml = [-70:1:70]/180*pi;
        angle_ml = [-180:-110]/180*pi;
        %angle_ml = 0:50;
        angle_sl1 = [-180:1:-180]/180*pi;
        angle_sl2 = [-110:1:180]/180*pi;
        angle_sl =[angle_sl1 angle_sl2];
        a_all_ref = a_all;
        p_ref = p;
        index_ml = find(scan_angle>=angle_ml(1) & scan_angle<=angle_ml(end));

        a_ml = a_all(:,index_ml);%main lobe steering vector
        p_ml = p(index_ml);      %main lobe expected response


        index_s1 = find(scan_angle>=angle_sl1(1) & scan_angle<=angle_sl1(end));     
        a_sl1  = a_all(:,index_s1);    %side lobe steering vector

        index_s2 = find(scan_angle>=angle_sl2(1) & scan_angle<=angle_sl2(end));    
        a_sl2  = a_all(:,index_s2);     %side lobe steering vector

        index_sl = [index_s1 index_s2];
        a_sl = [a_sl1 a_sl2];            %total sidelobe steering vector
        p_sl = p(index_sl);              %side lobe response
    end
    

    %figure(1)
    P_Hz=4000;
    P_B=round(P_Hz/6000*Nfft/2);
    %plot(scan_angle*180/pi,p,'b')
    %plot(scan_angle*180/pi,20*log10(p1(:,P_B)/max(p1(:,P_B))),'b')
    ylim([-50 0])
    xlabel('\phiAzimuth(^o)')
    ylabel('beam output/dB')
    title([num2str(P_Hz) ' Hzbeam pattern'])
    %p_all = [p_all p];
     figure(1)
    plot( scan_angle*180/pi, p,'-.','linewidth', 3)
    hold on; grid on;
    ylim([-50 0]);xlim([-180 180]);
    xlabel('angle/\circ','fontsize',14);
    ylabel('beam response/dB','fontsize',14);

    p_all = [p_all p];
    a0_all=[a0_all p11'];
    
end
legend({num2str(f_range(1)),num2str(f_range(2)),num2str(f_range(3)),...
        num2str(f_range(4)),num2str(f_range(5)),num2str(f_range(6)),...
        num2str(f_range(7))},'fontsize',14,'Location','best')
figure(5)
[ff,sacn_angle]=meshgrid(f_range(2:end),scan_angle*180/pi);
p_all(p_all<-90) = -90;
surf(sacn_angle,ff,p_all);
set(gca,'view',[135 45])
zlim([-90 0]) 
xlim([-180 180])
xlabel('angle/\circ','fontsize',14);
ylabel('frequency/Hz','fontsize',14);
zlabel('beam response/dB','fontsize',14);



 
 


p_all_opt=[];
Wopt_all=[];
w_msl_all=[];
for nf =2:Nfft/2+1
         freq = f_range(nf);
         %freq1=500;
    a0=exp(-1j*2*pi*R*freq/c*sin(seta0)*cos(fai0-gamma));   %计算期望导向矢量 %Calculate desired steering vector
    C=[a0];
    Rxx=squeeze(Fvv2(nf,:,:));
%   Rxx=1;
    Wc=(Rxx*C)/(C'*Rxx*C)*f;                % MVDR算法获得权系数
%       Zero point adjustment (support multiple zero points)
    Wopt=Wc; 
    step = 1;
   % scan_angle = [-90:step:90];                            %扫描角度范围 %Scan Angle Range    
    scan_angle=[-180:1:180]/180*pi; 
    a_all =[];                           
    for i_1 = 1:length(scan_angle) 
        as=exp(-1j*2*pi*R*freq/c*sin(seta)*cos(scan_angle(i_1)-gamma')); 
        p2(i_1,nf)=abs(as*conj(Wopt));
        a_all= [a_all as'];
    end
p = 20*log10(abs(a_all'*conj(Wopt))/max(abs(a_all'*conj(Wopt))));        %计算方向图 %Computational pattern
p2=Wopt';
%p = 20*log10(abs(a_all'*conj(a0))/max(abs(a_all'*conj(a0))));  
    %Side lobe control main lobe minimum error
    max_sl = -10; %%%
    max_sl=10.^(max_sl/20);
    p_ml = 10.^(p_ml/20);                           %%%期望主瓣响应 %expected main lobe response
    cvx_begin quiet
        variable w_msl(M) complex
        p_temp=(w_msl' * a_all);             %所有的乘积      % all the products       
        p_temp_ml=p_temp(index_ml);   %选出对应的主瓣 %select the coressponding main lobe
        p_temp_ml=p_temp_ml(:); 
        p_ml=p_ml(:);                              %期望主瓣响应 %expected main lobe response
        p_temp_sl=p_temp(index_sl);       %选出对应的旁瓣 %Select the corresponding side lobe

        minimize(norm(p_temp_ml-p_ml,2))  %最小化2范数  %Minimize the 2-norm
    subject to
        abs(p_temp_sl) <= max_sl;
        w_msl'*conj(a0) == 1;
        norm(w_msl,2) <=M;
    cvx_end
    p_opt = 20*log10(abs(w_msl' * a_all)/max(abs(w_msl' * a_all)));
    p111=w_msl';
    C1=[w_msl];
    Rxx1=squeeze(Fvv2(nf,:,:));
%   Rxx=1;
    Wc1=(Rxx1*C1)/(C1'*Rxx1*C1)*f;                % MVDR algorithm to obtain weight coefficient
%       Zero point adjustment (support multiple zero points)
    Wopt1=Wc1; 
    p_opt1 = 20*log10(abs((Wopt1)' * a_all)/max(abs((Wopt1)' * a_all)));

    %% 绘图 %drawing
    %figure(3)
    P_Hz=4000;
    P_B=round(P_Hz/6000*Nfft/2);
    %plot(scan_angle*180/pi, p_opt,'b')
    %plot(scan_angle*180/pi,20*log10(p1(:,P_B)/max(p1(:,P_B))),'b')
    ylim([-50 0])
    xlabel('\phiAzimuth(^o)')
    ylabel('beam output/dB')
    title([num2str(P_Hz) ' Hzbeam pattern'])
    p_all_opt = [p_all_opt p_opt.']; 
    w_msl_all=[w_msl_all p111'];
    Wopt_all=[Wopt_all p2'];
      figure(3)
    plot(scan_angle*180/pi,  p_opt,'-.','linewidth', 3)
    hold on; grid on;
    ylim([-50 0]);xlim([-180 180]);
    xlabel('angle/\circ','fontsize',14);
    ylabel('beam response/dB','fontsize',14);
end
legend({num2str(f_range(1)),num2str(f_range(2)),num2str(f_range(3)),...
        num2str(f_range(4)),num2str(f_range(5)),num2str(f_range(6)),...
        num2str(f_range(7))},'fontsize',14,'Location','best')
figure(4)
[sacn_angle,ff]=meshgrid(f_range(2:end),scan_angle*180/pi);
p_all_opt(p_all_opt<-180) = -180;
surf(ff,sacn_angle,p_all_opt);
set(gca,'view',[135 45])
zlim([-60 0]) 
xlim([-180 180])
xlabel('angle/^o','fontsize',14)
ylabel('frequency/Hz','fontsize',14)
zlabel('beam response/dB','fontsize',14)