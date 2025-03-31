function [e,w_hat,update_ratio] = Volterra_LMS1(x,d,imp,vi,Nw ,mu1)


%滤波参数设置
% Nw     = 128;  %系数个数
% mu=0.01;

dim=length(d);
w=zeros(Nw,dim);           % initial coefficient vector
w_hat=zeros(Nw,dim); 

y=zeros(dim,1);
e=zeros(dim,1);

time=0;
x_vec = zeros(Nw,1);
for i=1:dim
    x_vec = [x(i);x_vec(1:end-1)];
      y(i) = x_vec.'*w(:,i);
      e(i) = d(i) - y(i);
      w(:,i+1)=w(:,i)+2*mu1*e(i)*x_vec;

      time=time+1;
      w_hat(:,i)=w(:,i+1);
end
%emse
e=e-vi-imp; %列向量
update_ratio = time/dim;

