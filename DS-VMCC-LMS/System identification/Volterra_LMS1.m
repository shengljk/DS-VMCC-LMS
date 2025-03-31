function [e,w_hat,update_ratio] = Volterra_LMS1(x,d,imp,vi)


%滤波参数设置
Nw     = 9;  %系数个数
u      = diag([0.078*ones(1,3) 0.009*ones(1,Nw-3)]);  %收敛因子矩阵
% u = diag(0.0135*ones(1,Nw));

dim=length(d);
w=zeros(Nw,dim);           % initial coefficient vector
w_hat=zeros(Nw,dim); 

y=zeros(dim,1);
e=zeros(dim,1);

time=0;

for i=1:dim
      y(i)=w(:,i)'*x(:,i);  % output sample
      e(i)=d(i)-y(i);    % error sample  
                  
      w(:,i+1)=w(:,i)+2*u*e(i)*x(:,i);  % new coefficient vector
      time=time+1;

      w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp; %列向量
update_ratio = time/dim;

