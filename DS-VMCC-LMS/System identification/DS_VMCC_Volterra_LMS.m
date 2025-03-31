function [e,w_hat,update_ratio,detection,false_alarm] = DS_VMCC_Volterra_LMS(x,d,P_up,var_noise,imp,vi)

% 滤波器参数
L     = 9;  %系数个数
u      = diag([0.08*ones(1,3) 0.012*ones(1,L-3)]);  %收敛因子矩阵


dim=length(d);
w=zeros(L,dim);           % initial coefficient vector
w_hat=zeros(L,dim); 
y=zeros(dim,1);
e=zeros(dim,1);

% 数据选择参数
v = 5;
% mu = 1/(v*(length(L))*var(x));
alpha = (P_up/v)/(1-P_up/v);
tau = (1+alpha)*(qfuncinv(P_up/2))^2; 

%可变核宽度重要参数
beita=0.9;
deta0=40; %核宽度上限
c=40;     %常量

%误差方差估计
Nw2=11;     %误差向量（平方）观察窗
lambda = 0.9;
% C1 = 1.483*(1+5/(Nw2-1));
power_av = zeros(dim,0);
f = zeros(Nw2,1);

%误差信号估计值的绝对值
Nw1=20;     %误差向量（绝对值）观察窗
e_hat_abs = zeros(dim,0);    %误差信号估计值的绝对值
G = zeros(Nw1,1);
deta=zeros(1,dim);  %可变核宽度
C1 = 1.483*(1+5/(Nw1-1));

%MCC
Jc=zeros(1,dim);
tH=zeros(1,dim);
tL=zeros(1,dim);

time = 0;
time_d = 0;
time_f = 0;

% DS_LMS_MCC
for i = 1:dim
    y(i)=w(:,i)'*x(:,i);  % output sample
    e(i)=d(i)-y(i);    % error sample 

   %计算误差信号方差估计值
    if i == 1
        power_av(i) = (e(i)^2);  %注意
    else
        f = [(e(i)^2); f(1:end-1)];
        power_av(i) = lambda*power_av(i-1) + (1-lambda)*median(f);
    end

  %计算误差信号估计值的绝对值
    if i==1
        e_hat_abs(i)=abs(e(i));
    else
        G=[abs(e(i));G(1:end-1)];
        e_hat_abs(i)=beita*abs(e_hat_abs(i-1))+C1*(1-beita)*min(G);
    end

    %计算可变核宽度的值
    if e_hat_abs(i) > (deta0/c)      %核宽度达到上限
        deta(i)=deta0;
    else
        deta(i)=c*e_hat_abs(i);
    end

    %计算MCC、tH、tL
    Jc(i)=exp(e(i)^2/(-2*deta(i)^2));
    tH(i)=exp(tau*var_noise/(-2*deta(i)^2));
%     tH(i)=exp(2*var_noise/(-2*deta(i)^2));
    tL(i)=exp((L+1)*power_av(i)/(-2*deta(i)^2));

    % Update Equation
    if Jc(i)>tH(i)    %创新不足
        w(:,i+1)=w(:,i);
    elseif Jc(i)>=tL(i)
        w(:,i+1)=w(:,i)+2*u*e(i)*x(:,i);
        time=time+1;   %更新次数
    else
        w(:,i+1)=w(:,i);
        if imp(i) ~= 0 %實際有:detection
                time_d = time_d + 1;   %正确检测脉冲次数
        else %實際無:false alarm
                time_f = time_f + 1;
        end
    end
     w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp;
update_ratio = time/dim;
detection = (time_d)/(sum(imp~=0));
false_alarm = (time_f)/(sum(imp==0));










