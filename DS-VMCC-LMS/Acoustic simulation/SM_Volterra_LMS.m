function [e,w_hat,update_ratio] = SM_Volterra_LMS(x,d,P_up,var_noise,imp,vi,L,mu2)


% 滤波器参数
% L     = 128;  %系数个数
% mu2=0.01;

dim=length(d);
w=zeros(L,dim);           % initial coefficient vector

w_hat=zeros(L,dim); 

y=zeros(dim,1);
e=zeros(dim,1);

time = 0;
x_vec = zeros(L,1);
% DS_LMS
for i=1:dim
    x_vec = [x(i);x_vec(1:end-1)];
    y(i) = x_vec.'*w(:,i);
    e(i) = d(i) - y(i);
    
    % Update Equation
    if abs(e(i))^2 <5*var_noise  %下限，上限为∞，误差小于下限阈值，输入创新不足
        w(:,i+1)=w(:,i);
%     if  abs(e(i))^2 <= 3*var_noise
%           w(:,i+1)=w(:,i);
    else        
            w(:,i+1)=w(:,i)+2*mu2*e(i)*x_vec;
            time = time + 1;  %权重更新次数加1
    end
     w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp;
update_ratio = time/dim;

