function [e,w_hat,update_ratio,detection,false_alarm] = DS_d1_Volterra_LMS(x,d,P_up,var_noise,imp,vi,L,mu3)

% % Calculate threshold
% v = 5;
% alpha = (P_up/v)/(1-P_up/v);
% tau = (1+alpha)*(qfuncinv(P_up/2))^2;

% 滤波器参数
% L     = 128;  %系数个数
% mu3=0.01;

dim=length(d);
w=zeros(L,dim);           % initial coefficient vector
w_hat=zeros(L,dim); 
y=zeros(dim,1);
e=zeros(dim,1);

% Average Power Estimation Parameter
Nw = 11;
lambda = 0.9;
C1 = 1.483*(1+5/(Nw-1));
power_av = zeros(dim,0);
f = zeros(Nw,1);

time = 0;
time_d = 0; %脉冲检测次数
time_f = 0; %脉冲检测错误次数

x_vec = zeros(L,1);
for i = 1:dim
   x_vec = [x(i);x_vec(1:end-1)];
   y(i) = x_vec.'*w(:,i);
   e(i) = d(i) - y(i);

    if i==1
        power_av(i) = (d(i)^2);
    else
        f = [(d(i)^2); f(1:end-1)];
        power_av(i) = lambda*power_av(i-1) + C1*(1-lambda)*median(f);
    end

   if (d(i)^2) > (L+1)*power_av(i) %判斷有imp
            w(:,i+1)=w(:,i);
            if imp(i) ~= 0 %實際有:detection
                time_d = time_d + 1;
            else %實際無:false alarm
                time_f = time_f + 1;
            end
   else %判斷沒有
            w(:,i+1)=w(:,i)+2*mu3*e(i)*x_vec;
            time = time + 1;  %权重更新次数加1
   end

    w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp;
update_ratio = time/dim;
detection = (time_d)/(sum(imp~=0));
false_alarm = (time_f)/(sum(imp==0));








