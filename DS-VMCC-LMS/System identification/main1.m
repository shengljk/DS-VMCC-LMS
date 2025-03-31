clear;	
%  高斯白输入
% Input: 
Nr     = 100;   %模拟次数
L=9;            %系数个数
dim    = 4e3;   %迭代次数，输入信号长度
var_x  = 1;    %输入信号方差


IN_MODE=2 ; %1 weak   2 strong

%核向量
w = [-0.76,-1.0,1.0,0.5,0,2.0,-1.6,0.8,1.2]';

EMSE1=zeros(Nr,dim);
EMSE_diniz=zeros(Nr,dim);
EMSE_d1=zeros(Nr,dim);
EMSE_MCC1=zeros(Nr,dim);
EMSE_MCC2=zeros(Nr,dim);

NMSD1= zeros(Nr,dim);
NMSD_diniz= zeros(Nr,dim);
NMSD_d1=zeros(Nr,dim);
NMSD_MCC1=zeros(Nr,dim);
NMSD_MCC2=zeros(Nr,dim);

W_1=zeros(L,Nr);
W_diniz=zeros(L,Nr);
W_d1=zeros(L,Nr);
W_MCC1=zeros(L,Nr);
W_MCC2=zeros(L,Nr);



w0_1=zeros(Nr,dim);
w0_diniz=zeros(Nr,dim);
w0_d1=zeros(Nr,dim);
w0_MCC1=zeros(Nr,dim);
w0_MCC2=zeros(Nr,dim);

w00_1=zeros(Nr,dim);
w00_diniz=zeros(Nr,dim);
w00_d1=zeros(Nr,dim);
w00_MCC1=zeros(Nr,dim);
w00_MCC2=zeros(Nr,dim);

w_1=zeros(Nr,dim);
w_diniz=zeros(Nr,dim);
w_d1=zeros(Nr,dim);
w_MCC1=zeros(Nr,dim);
w_MCC2=zeros(Nr,dim);

for j=1:Nr

%    （1）获得输入信号
   x=sqrt(var_x)*randn(dim,1);  %高斯白输入（均值为0，方差为1）
  
   xl1=zeros(dim,1); xl2=xl1; 
   xl1(2:dim)=x(1:dim-1);    % x(k-1)
   xl2(3:dim)=x(1:dim-2);    % x(k-2)
   uxl=[x xl1 xl2 x.^2 x.*xl1 x.*xl2 xl1.^2 xl1.*xl2 xl2.^2]'; % 输入信号最终形式

   %（2）获得期望信号
   y=-.76*x-xl1+xl2+.5*x.^2+2*x.*xl2-1.6*xl1.^2+1.2*xl2.^2+.8*xl1.*xl2;  %非线性系统输出

   % 测量噪声
%     var_y = var(y);             %未知系统输出方差
%     SNRdB = 20;
%     var_v = var_y./(10^(SNRdB/10));   
    var_v=0.01;
    vi = sqrt(var_v).*randn(dim,1);   %加性高斯白噪声输入
    dn1 = y + vi;           % dn1期望信号（只有加性高斯白噪声）

    % Impulse noise
    switch(IN_MODE)      
        case 1 %weak
            disp('Weak IN');
            GINR = 0.01;
            pb = 0.05;
        case 2 %strong
            disp('Strongly IN');
            GINR = 0.001;
            pb = 0.1;
        otherwise
            disp('ERROR IN MODE');
            exit();
    end
    disp(j)

    sigma=sqrt(var_v);
    imp = BG_Noise(pb, sigma ,GINR,length(y));   %产生脉冲噪声BG建模
    dn2 = dn1 + imp;     % dn2期望信号（加性高斯白噪声和脉冲噪声）  


    %（3）代入滤波算法
    P_up=0.3;

    [e,w1_hat]=Volterra_LMS1(uxl,dn2,imp,vi);
    [e_diniz,w_diniz_hat,update_ratio1(j,1)] = SM_Volterra_LMS(uxl,dn2,P_up,var_v,imp,vi);
    [e_MCC1,w_MCC1_hat,update_ratio3(j,1),detection2(j,1),false_alarm2(j,1)] = MCC_Volterra_LMS(uxl,dn2,P_up,var_v,imp,vi);
    [e_MCC2,w_MCC2_hat,update_ratio4(j,1),detection3(j,1),false_alarm3(j,1)] = DS_VMCC_Volterra_LMS(uxl,dn2,P_up,var_v,imp,vi);
    [e_d1,w_d1_hat,update_ratio5(j,1),detection4(j,1),false_alarm4(j,1)] = DS_d1_Volterra_LMS(uxl,dn2,P_up,var_v,imp,vi);


    %均方误差
    EMSE1(j,:)=e'.^2;
    EMSE_diniz(j,:)=e_diniz'.^2;
    EMSE_d1(j,:)=e_d1'.^2;
    EMSE_MCC1(j,:)=e_MCC1'.^2;
    EMSE_MCC2(j,:)=e_MCC2'.^2;

%系数
W_1(:,j)=w1_hat(:,end);
W_diniz(:,j)=w_diniz_hat(:,end);
W_d1(:,j)=w_d1_hat(:,end);
W_MCC1(:,j)=w_MCC1_hat(:,end);
W_MCC2(:,j)=w_MCC2_hat(:,end);

%
NMSD1(j,:) = Normalized_Mean_Square_Deviation2(w,w1_hat);
NMSD_diniz(j,:)=Normalized_Mean_Square_Deviation2(w,w_diniz_hat);
NMSD_d1(j,:)=Normalized_Mean_Square_Deviation2(w,w_d1_hat);
NMSD_MCC1(j,:)=Normalized_Mean_Square_Deviation2(w,w_MCC1_hat);
NMSD_MCC2(j,:)=Normalized_Mean_Square_Deviation2(w,w_MCC2_hat);

%权系数 w0(k)
w0_1=w1_hat(1,:);
w0_diniz=w_diniz_hat(1,:);
w0_d1=w_d1_hat(1,:);
w0_MCC1=w_MCC1_hat(1,:);
w0_MCC2=w_MCC2_hat(1,:);

% w0,0(k)
w00_1=w1_hat(4,:);
w00_diniz=w_diniz_hat(4,:);
w00_d1=w_d1_hat(4,:);
w00_MCC1=w_MCC1_hat(4,:);
w00_MCC2=w_MCC2_hat(4,:);

% 其他权值，如改变w1_hat的行，可得不同权重
w_1=w1_hat(3,:);
w_diniz=w_diniz_hat(3,:);
w_d1=w_d1_hat(3,:);
w_MCC1=w_MCC1_hat(3,:);
w_MCC2=w_MCC2_hat(3,:);

% 噪声方差
% VN(j)=(1-pb)*(1+update_ratio4(j,1)/(5-update_ratio4(j,1)))*var_v+pb*(1-detection3(j,1))*var(e_MCC2);
VN(j)=(1-pb)*var_v+pb*(1-detection3(j,1))*var(e_MCC2);

end

% 噪声方差
VN_av=mean(VN);

EMSE1_av=10*log10(sum(EMSE1,1)/Nr);  
EMSE_diniz_av=10*log10(sum(EMSE_diniz,1)/Nr);
EMSE_d1_av=10*log10(sum(EMSE_d1,1)/Nr); 
EMSE_MCC1_av=10*log10(sum(EMSE_MCC1,1)/Nr);
EMSE_MCC2_av=10*log10(sum(EMSE_MCC2,1)/Nr);

NMSD1 = 10*log10(sum(NMSD1,1)/Nr);
NMSD_diniz=10*log10(sum(NMSD_diniz,1)/Nr);
NMSD_d1=10*log10(sum(NMSD_d1,1)/Nr);
NMSD_MCC1=10*log10(sum(NMSD_MCC1,1)/Nr);
NMSD_MCC2=10*log10(sum(NMSD_MCC2,1)/Nr);


% 求辨识参数
W_1_av=sum(W_1,2)/Nr;
W_diniz_av=sum(W_diniz,2)/Nr;
W_d1_av=sum(W_d1,2)/Nr;
W_MCC1_av=sum(W_MCC1,2)/Nr;
W_MCC2_av=sum(W_MCC2,2)/Nr;

% 核参数w0(k)
w0_1_av=mean(w0_1,1);
w0_diniz_av=mean(w0_diniz,1);
w0_d1_av=mean(w0_d1,1);
w0_MCC1_av=mean(w0_MCC1,1);
w0_MCC2_av=mean(w0_MCC2,1);

w00_1_av=mean(w00_1,1);
w00_diniz_av=mean(w00_diniz,1);
w00_d1_av=mean(w00_d1,1);
w00_MCC1_av=mean(w00_MCC1,1);
w00_MCC2_av=mean(w00_MCC2,1);

w_1_av=mean(w_1,1);          
w_diniz_av=mean(w_diniz,1);
w_d1_av=mean(w_d1,1);
w_MCC1_av=mean(w_MCC1,1);
w_MCC2_av=mean(w_MCC2,1);

% 理论MSD
R = (uxl * uxl') / dim; % 输入信号的协方差矩阵
u      = diag([0.08*ones(1,3) 0.012*ones(1,L-3)]); 
tr_u = trace(u); % tr(u)
tr_uR = trace(u * R); % tr(u * R)
p_up=mean(update_ratio4);
% 计算稳态 MSD
MSD = (tr_u *p_up*VN_av*L)/(2-p_up*tr_uR);
NMSD_theory = 10 * log10(MSD / norm(w)^2); % 归一化MSD



% 辨识参数
disp('***************************************');
disp(['（VLMS） :   ',num2str(W_1_av')]);
disp(['（SM-VLMS） :   ',num2str(W_diniz_av')]);
disp(['（DS_Jenong_VLMS） :   ',num2str(W_d1_av')]);  
disp(['（MCCDS_VLMS） :   ',num2str(W_MCC1_av')]);
disp(['（Proposed） :   ',num2str(W_MCC2_av')]);

% %稳态MSE
disp('***************************************');
disp(['（VLMS） :   ',num2str(EMSE1_av(end))]);
disp(['（SM-VLMS） :   ',num2str(EMSE_diniz_av(end))]);
disp(['（DS_Jenong_VLMS） :   ',num2str(EMSE_d1_av(end))]);
disp(['（MCCDS_VLMS） :   ',num2str(EMSE_MCC1_av(end))]);
disp(['（Proposed） :   ',num2str(EMSE_MCC2_av(end))]);

% 稳态NMSD
disp('***************************************');
disp(['（VLMS） :   ',num2str(NMSD1(end))]);
disp(['（SM-VLMS） :   ',num2str(NMSD_diniz(end))]);
disp(['（DS_Jenong_VLMS） :   ',num2str(NMSD_d1(end))]);
disp(['（MCCDS_VLMS） :   ',num2str(NMSD_MCC1(end))]);
disp(['（Proposed） :   ',num2str(NMSD_MCC2(end))]);
disp(['（Proposed-Theoretical） :   ',num2str(NMSD_theory)]);



% % %更新率
disp('***************************************');
disp(['（SM-VLMS） Update Ratio:   ',num2str(mean(update_ratio1))]);
disp(['（DS_Jenong_VLMS） Update Ratio:   ',num2str(mean(update_ratio5))]);
disp(['（MCCDS_VLMS） Update Ratio:   ',num2str(mean(update_ratio3))]);
disp(['（Proposed） Update Ratio:   ',num2str(mean(update_ratio4))]);

% 检测率
disp('***************************************');
disp(['（DS_Jenong_Volterra_LMS）  Detection Rate:   ',num2str(mean(detection4))]);
disp(['（MCCDS_Volterra_LMS）  Detection Rate:   ',num2str(mean(detection2))]);
disp(['（Proposed）  Detection Rate:   ',num2str(mean(detection3))]);

%误判率
disp('***************************************');
disp(['（DS_Jenong_Volterra_LMS）  false alarm rate:   ',num2str(mean(false_alarm4))]);
disp(['（MCCDS_Volterra_LMS）  false alarm rate:   ',num2str(mean(false_alarm2))]);
disp(['（Proposed）  false alarm rate:   ',num2str(mean(false_alarm3))]);


% 绘图

light_colors = [
    1, 0, 0;   % 红色
    0, 1, 0;   % 绿色
    0, 0, 1;   % 蓝色
    1, 1, 0;   % 黄色
    1, 0, 1;   % 品红色
    0, 1, 1;   % 青色
    1, 0.5, 0; % 橙色
];

figure,
plot(1:dim,EMSE1_av, 'color',light_colors(2, :), 'linewidth', 1);
hold on
plot(1:dim,EMSE_diniz_av,'color',light_colors(6, :),'LineWidth',1);
plot(1:dim,EMSE_d1_av,'color',light_colors(3, :),'LineWidth',1);
plot(1:dim,EMSE_MCC1_av,'color',light_colors(4, :),'LineWidth',1);
plot(1:dim,EMSE_MCC2_av,'color',light_colors(1, :),'LineWidth',1);

title('Learning Curve for EMSE');
xlabel('Number of iterations, k'); ylabel('EMSE [dB]');
legend('VLMS','SM-VLMS','DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS');

figure,
plot(1:dim,NMSD1, 'color',light_colors(2, :),'linewidth', 1);
hold on
plot(1:dim,NMSD_diniz,'color',light_colors(6, :),'LineWidth',1);
plot(1:dim,NMSD_d1,'color',light_colors(3, :),'LineWidth',1);
plot(1:dim,NMSD_MCC1,'color',light_colors(4, :),'LineWidth',1);
plot(1:dim,NMSD_MCC2,'color',light_colors(1, :),'LineWidth',1);
plot(1:dim,NMSD_theory*ones(1,dim),'--k');

title('Learning Curve for NMSD');
xlabel('Number of iterations, k'); ylabel('NMSD [dB]');
legend('VLMS','SM-VLMS','DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS','Theoretical value');


%% 权值收敛曲线  可以改成2000次迭代
% figure,
% plot(1:dim,w0_d1_av,'color',light_colors(3, :),'LineWidth',1);
% hold on
% plot(1:dim,w0_MCC1_av,'color',light_colors(4, :),'LineWidth',1);
% plot(1:dim,w0_MCC2_av,'color',light_colors(1, :),'LineWidth',1);
% plot(1:dim,-0.76*ones(1,dim),'--k',"LineWidth",1.5);
% title('Weight Convergence Curve');
% xlabel('Number of iterations, k'); ylabel('w_0(k)');
% legend('DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS','h_0(k)');
% axis([0 dim -2 2]);


% figure,
% plot(1:dim,w00_d1_av,'color',light_colors(3, :),'LineWidth',1);
% hold on
% plot(1:dim,w00_MCC1_av,'color',light_colors(4, :),'LineWidth',1);
% plot(1:dim,w00_MCC2_av,'color',light_colors(1, :),'LineWidth',1);
% plot(1:dim,0.5*ones(1,dim),'--k',"LineWidth",1.5);
%  title('Weight Convergence Curve');
% xlabel('Number of iterations, k'); ylabel('w_{0,0}(k)');
% legend('DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS','h_{0,0}(k)');  
% axis([0 dim -0.5 2]);

% figure,                                                              %其他权值
% plot(1:dim,w_d1_av,'color',light_colors(3, :),'LineWidth',1);
% hold on
% plot(1:dim,w_MCC1_av,'color',light_colors(4, :),'LineWidth',1);
% plot(1:dim,w_MCC2_av,'color',light_colors(1, :),'LineWidth',1);
% plot(1:dim,1*ones(1,dim),'--k',"LineWidth",1.5);
%  title('Weight Convergence Curve');
% xlabel('Number of iterations, k'); ylabel('w_{2}(k)');
% legend('DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS','h_{2}(k)');  
% axis([0 dim 0 3]);


