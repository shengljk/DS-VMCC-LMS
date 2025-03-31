function [e,w_hat,update_ratio] = SM_Volterra_LMS(x,d,P_up,var_noise,imp,vi)

%数据选择参数
v=5;
alpha = (P_up/v)/(1-P_up/v);
tau = (1+alpha)*(qfuncinv(P_up/2))^2;

% 滤波器参数
L     = 9;  %系数个数
u      = diag([0.078*ones(1,3) 0.009*ones(1,L-3)]);  %收敛因子矩阵

dim=length(d);
w=zeros(L,dim);           % initial coefficient vector

w_hat=zeros(L,dim); 

y=zeros(dim,1);
e=zeros(dim,1);

time = 0;

% DS_LMS
for i=1:dim
    y(i)=w(:,i)'*x(:,i);  % output sample
    e(i)=d(i)-y(i);    % error sample 
    
    % Update Equation
    if abs(e(i))^2 <5*var_noise  %下限，上限为∞，误差小于下限阈值，输入创新不足
        w(:,i+1)=w(:,i);
%     if  abs(e(i))^2 <= 3*var_noise
%           w(:,i+1)=w(:,i);
    else        
            w(:,i+1)=w(:,i)+2*u*e(i)*x(:,i);
            time = time + 1;  %权重更新次数加1
    end
     w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp;
update_ratio = time/dim;

