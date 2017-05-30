function [R]=Deltak(Tau, epsilon0, delta0, sigma);
% 
% calculate RSN * Min(1/(sqrt(1/N_k)+(1/N_k+1))
%  whis provides an estimation of the shift of mean value of t-values 
% of true positive with respect ot the false alarms

for k=2:length(Tau)
    N(k)=Tau(k)-Tau(k-1)-2*epsilon0;
end;

for k=1:length(Tau)-2
    R(k)=(delta0/sigma)/sqrt(1/N(k+1)+1/N(k+2));
end;

end

