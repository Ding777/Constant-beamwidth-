clc
close all;
warning off;
clear all; close all; clc;
%% Basic parameters
f01 = 4000;           % Signal center frequency
B =  8000; 
Nfft =512;                                % FFT points
F0=linspace(0, 8000, Nfft/2+1);           % Actual frequency corresponding to each FFT bin
f_range=linspace(f01-B/2, f01+B/2, Nfft/2+1);     % Frequency range around center frequency
% F(1:10)=F(11);
% F(1:3)=F(4);

% F(1:16)=F(17);
% F(129:257)=F(128);
% F(1:end)=F(70);

% F(200:257)=F(199);
% F(1:128)=F(129);
% F = fl:B/(Nfft/2):fh; 

M =16;                                % Number of array elements
c = 340;                              % Speed of sound (m/s)
lambda = c/f01;                        % Center wavelength
R =0.025;                               % Actual radius of the circular array
% R =0.5*M*lambda/(2*pi);             % Ideal radius (example)

gamma=2*pi*[0:M-1]'/M ;
% d = 0.5*lambda;                        % Element spacing (half-wavelength ideal)
% BW= 0.886*lambda/(M*d)*180/pi;         % Half-power beamwidth at center frequency (degrees)
% fprintf('Half-power beamwidth BW = %0.2f°\n',BW)

% seta1=[60]*pi/180;                     % Zero-constraint elevation angle 1
% fai1=[-140]*pi/180;                    % Zero-constraint azimuth 1  
 
% seta2=[60]*pi/180;                     % Zero-constraint elevation angle 2
% fai2=[-111]*pi/180;                    % Zero-constraint azimuth 2 

%% Read actual audio / parameters for DOA
seta0 = 60*pi/180;                       % Signal arrival elevation angle (radians)
fai0  =-150*pi/180;                      % Signal arrival azimuth angle (radians)
f=[1]';    
%% Diffuse noise field spatial coherence (MSC)
Fs=16000;                                                              
f0 = (1:Nfft/2+1)*Fs/Nfft;
Fvv = zeros(Nfft/2+1,M,M);
 %f0=500;
% dij = 2*R*sin(pi/M);

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

% Diagonal loading mapping curve (for stability)
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
%M = 10;                 % Number of elements
%dd= (0:1:M-1) * (c/f0/2);   % Element positions
p_all=[];
a0_all=[];
for nf =2:Nfft/2+1
    freq = f_range(nf);
    freq1=500;
   % a0 = exp(-1i * 2 * pi * freq /c *dd.'*sind(theta0));   % Compute desired steering vector
    a0=exp(-1j*2*pi*R*freq/c*sin(seta0)*cos(fai0-gamma));   % Calculate desired steering vector (UCA)
    a0=a0;
    C=[a0];
   % a0=a0/100;
    Rxx=squeeze(Fvv2(nf,:,:));
%   Rxx=1;
    Wc=(Rxx*C)/(C'*Rxx*C)*f;                % MVDR algorithm to obtain weight coefficients
%       Null steering (supports multiple nulls)
    Wopt=Wc;  
    step = 1;
    seta=seta0; 
    scan_angle = [-180:step:180];                        % Scan angle range (degrees)
    scan_angle=[-180:step:180]/180*pi;                   % convert to radians
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
    p = 20*log10(abs(a_all'*conj(Wopt))/max(abs(a_all'*conj(Wopt))));  % Compute beam pattern (dB)
    p3 = 20*log10(abs(a_all')/max(abs(a_all')));  
    p4 = 20*log10(abs(p1)/max(abs(p1))); 
    %p = 20*log10(abs(a_all'*conj(a0))/max(abs(a_all'*conj(a0))));  
    p7 = 20*log10(abs(a_all'*(Wopt))/max(abs(a_all'*(Wopt))));
    p11=conj(a0)';
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

        a_ml = a_all(:,index_ml); % Main-lobe steering vectors
        p_ml = p(index_ml);       % Expected main-lobe response


        index_s1 = find(scan_angle>=angle_sl1(1) & scan_angle<=angle_sl1(end));     
        a_sl1  = a_all(:,index_s1);    % Side-lobe steering vectors

        index_s2 = find(scan_angle>=angle_sl2(1) & scan_angle<=angle_sl2(end));    
        a_sl2  = a_all(:,index_s2);    % Side-lobe steering vectors

        index_sl = [index_s1 index_s2];
        a_sl = [a_sl1 a_sl2];          % Total sidelobe steering vectors
        p_sl = p(index_sl);            % Side-lobe responses
    end
    

    %figure(1)
    P_Hz=4000;
    P_B=round(P_Hz/6000*Nfft/2);
    %plot(scan_angle*180/pi,p,'b')
    %plot(scan_angle*180/pi,20*log10(p1(:,P_B)/max(p1(:,P_B))),'b')
    ylim([-50 0])
    xlabel('\phiAzimuth(^o)')
    ylabel('beam output/dB')
    title([num2str(P_Hz) ' Hz beam pattern'])
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
    a0=exp(-1j*2*pi*R*freq/c*sin(seta0)*cos(fai0-gamma));   % Calculate desired steering vector
    C=[a0];
    Rxx=squeeze(Fvv2(nf,:,:));
%   Rxx=1;
    Wc=(Rxx*C)/(C'*Rxx*C)*f;                % MVDR algorithm to obtain weight coefficients
%       Zero point adjustment (support multiple nulls)
    Wopt=Wc; 
    step = 1;
   % scan_angle = [-90:step:90];                            % Scan angle range    
    scan_angle=[-180:1:180]/180*pi; 
    a_all =[];                           
    for i_1 = 1:length(scan_angle) 
        as=exp(-1j*2*pi*R*freq/c*sin(seta)*cos(scan_angle(i_1)-gamma')); 
        p2(i_1,nf)=abs(as*conj(Wopt));
        a_all= [a_all as'];
    end
p = 20*log10(abs(a_all'*conj(Wopt))/max(abs(a_all'*conj(Wopt))));        % Compute beam pattern (dB)
p2=Wopt';
%p = 20*log10(abs(a_all'*conj(a0))/max(abs(a_all'*conj(a0))));  
    % Sidelobe control with main-lobe least-squares
    max_sl = -10; % dB threshold for sidelobe
    max_sl=10.^(max_sl/20);
    p_ml = 10.^(p_ml/20);                           % Expected main-lobe response (linear)
    cvx_begin quiet
        variable w_msl(M) complex
        p_temp=(w_msl' * a_all);             % All products       
        p_temp_ml=p_temp(index_ml);   % Select corresponding main-lobe entries
        p_temp_ml=p_temp_ml(:); 
        p_ml=p_ml(:);                       % Expected main-lobe response (vector)
        p_temp_sl=p_temp(index_sl);         % Select corresponding sidelobe entries

        minimize(norm(p_temp_ml-p_ml,2))  % Minimize 2-norm (least squares)
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
    Wc1=(Rxx1*C1)/(C1'*Rxx1*C1)*f;                % MVDR algorithm to obtain weight coefficients
%       Zero point adjustment (support multiple nulls)
    Wopt1=Wc1; 
    p_opt1 = 20*log10(abs((Wopt1)' * a_all)/max(abs((Wopt1)' * a_all)));

    %% Plotting
    %figure(3)
    P_Hz=4000;
    P_B=round(P_Hz/6000*Nfft/2);
    %plot(scan_angle*180/pi, p_opt,'b')
    %plot(scan_angle*180/pi,20*log10(p1(:,P_B)/max(p1(:,P_B))),'b')
    ylim([-50 0])
    xlabel('\phiAzimuth(^o)')
    ylabel('beam output/dB')
    title([num2str(P_Hz) ' Hz beam pattern'])
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
