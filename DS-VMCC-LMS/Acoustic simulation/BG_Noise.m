function y = BG_Noise(p, sigma ,GINR,N)

% p : probability of impulsive noise occurrence 
% sigma  : standardized deviation of AWGN
% GINR : Gaussian to Impulsive noise ratio
% N : length of noise

b = binornd(1,p,N,1);
g = randn(N,1)*sigma*sqrt(( 1/GINR));      %加上根号
y = (b.*g);  %噪声
