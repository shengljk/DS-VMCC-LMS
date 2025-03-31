clear;clc;

% 读取音频文件
filename = 'clean_signal_audio2.wav'; 
[y, fs] = audioread(filename);

% 计算样本范围
start_time = 5; % 起始时间 (秒)
end_time = 12; % 结束时间 (秒)
start_sample = round(start_time * fs) + 1; % 起始样本索引
end_sample = round(end_time * fs); % 结束样本索引

% 读取部分片段
y_segment = y(start_sample:end_sample, :); % 读取指定范围的样本(第5到7s)
% 输入语音
x=y_segment;

% Input: 
Nr     = 10; % 模拟次数
dim    = length(x);

var_v  =0.01;

% % 声学通道脉冲响应
N =128; % 脉冲响应长度
alpha = 0.9; % 衰减因子
% 生成指数衰减脉冲响应
w0 = alpha.^(0:N-1);
w0=w0';

% 步长
mu1     = 0.005;
mu2     = 0.005;
mu3     = 0.005;
mu4     = 0.005;
mu5     = 0.005;

NMSD1= zeros(Nr,dim);
NMSD_diniz= zeros(Nr,dim);
NMSD_d1=zeros(Nr,dim);
NMSD_MCC1=zeros(Nr,dim);
NMSD_MCC2=zeros(Nr,dim);

% Body:
for j=1:Nr
   disp(j) 
   vi = sqrt(var_v).*randn(dim,1);   %加性高斯白噪声输入        
   x=x';                     % 输入
   w=zeros(N,dim);           % initial coefficient vector
   y=zeros(dim,1);
   % Impulse noise
    GINR = 0.0005;
    pb = 0.1;
    sigma=sqrt(var_v);
    imp = BG_Noise(pb, sigma ,GINR,length(y));   %产生脉冲噪声BG建模
   % 求期望
   y1 = zeros(dim,1);   %未知系统输出
   x_vec1 = zeros(N,1);
    for i = 1:dim
        x_vec1 = [x(i); x_vec1(1:end-1)];
        y1(i) = x_vec1.'*w0;   
    end
    d=y1+vi+imp;
  

     %（3）代入滤波算法
    P_up=0.4;

    [e,w1_hat]=Volterra_LMS1(x,d,imp,vi,N,mu1);
    [e_diniz,w_diniz_hat,update_ratio1(j,1)] = SM_Volterra_LMS(x,d,P_up,var_v,imp,vi,N,mu2);
    [e_MCC1,w_MCC1_hat,update_ratio3(j,1),detection2(j,1),false_alarm2(j,1)] = MCC_Volterra_LMS(x,d,P_up,var_v,imp,vi,N,mu3);
    [e_MCC2,w_MCC2_hat,update_ratio4(j,1),detection3(j,1),false_alarm3(j,1),Y] = DS_VMCC_Volterra_LMS(x,d,P_up,var_v,imp,vi,N,mu4);
    [e_d1,w_d1_hat,update_ratio5(j,1),detection4(j,1),false_alarm4(j,1)] = DS_d1_Volterra_LMS(x,d,P_up,var_v,imp,vi,N,mu5);

    %
NMSD1(j,:) = Normalized_Mean_Square_Deviation2(w0,w1_hat);
NMSD_diniz(j,:)=Normalized_Mean_Square_Deviation2(w0,w_diniz_hat);
NMSD_d1(j,:)=Normalized_Mean_Square_Deviation2(w0,w_d1_hat);
NMSD_MCC1(j,:)=Normalized_Mean_Square_Deviation2(w0,w_MCC1_hat);
NMSD_MCC2(j,:)=Normalized_Mean_Square_Deviation2(w0,w_MCC2_hat);  

% 还原后语音
YY(j,:)=Y;

end

NMSD1 = 10*log10(sum(NMSD1,1)/Nr);
NMSD_diniz=10*log10(sum(NMSD_diniz,1)/Nr);
NMSD_d1=10*log10(sum(NMSD_d1,1)/Nr);
NMSD_MCC1=10*log10(sum(NMSD_MCC1,1)/Nr);
NMSD_MCC2=10*log10(sum(NMSD_MCC2,1)/Nr);

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
plot(1:dim,NMSD1, 'color',light_colors(2, :),'linewidth', 1);
hold on
plot(1:dim,NMSD_diniz,'color',light_colors(6, :),'LineWidth',1);
plot(1:dim,NMSD_d1,'color',light_colors(3, :),'LineWidth',1);
plot(1:dim,NMSD_MCC1,'color',light_colors(4, :),'LineWidth',1);
plot(1:dim,NMSD_MCC2,'color',light_colors(1, :),'LineWidth',1);
title('Learning Curve for NMSD');
xlabel('Number of iterations, k'); ylabel('NMSD [dB]');
legend('VLMS','SM-VLMS','DS-Jeong-VLMS','MCCDS-VLMS','DS-VMCC-VLMS');

% 播放语音信号
% sound(y, fs);
% figure,
% t = (0:length(x)-1)/fs;
% plot(t, x);
% xlabel('时间 (秒)');
% ylabel('幅度');
% title('加载的语音信号');

figure;

% 绘制输入语音信号
subplot(3,1,1);
plot(x,'-k');
title('Input Voice Signal');


% 绘制带有噪声的回波信号
subplot(3,1,2);
plot(d,'-k');
title('Echo Voice Signal (with Noise)');


% 绘制恢复的语音信号
subplot(3,1,3);
plot(mean(YY),'-k');
title('Restored Voice Signal');
