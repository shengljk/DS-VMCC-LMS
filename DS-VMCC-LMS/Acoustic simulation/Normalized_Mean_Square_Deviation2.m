function NMSD = Normalized_Mean_Square_Deviation2(w,w_hat)

N = size(w_hat,2);
L = length(w);

for i = 1:N
    
        NMSD(i) = norm(w-w_hat(:,i))^2/norm(w)^2;
      
end