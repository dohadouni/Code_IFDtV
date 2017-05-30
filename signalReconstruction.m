function [hat_s] = signalReconstruction(hat_Tau, hat_mu)
%% Description of the function:
%
% Reconstruct the signal hat_s, as a piecewise constant function
% with jumps at hat_tau and values hat_mu(k+1) between two successive jumps
%
% 
% INPUT
% hat_Tau: family of K change points, with last value hat_Tau(K)=length(X)
% and 'first value' hat_Tau(0)=0
% hat_mu: family of K mean value hat_mu(k) for
% i\in[hat_Tau(k)+1:hat_Tau(k+1)]. 
% 
% WARNING:
% we should have length(hat_mu)= length(hat_Tau)
%

% hat_mu = [hat_mu hat_mu(length(hat_mu))];
for i=1:hat_Tau(1)
    hat_s(i)=hat_mu(1);
end;



for k=1:length(hat_Tau)-1
    for i=hat_Tau(k)+1:hat_Tau(k+1)
        hat_s(i)=hat_mu(k+1);
    end;
end;
% 
% for i=hat_Tau(length(hat_Tau)-1):hat_Tau(length(hat_Tau))
%     hat_s(i) = hat_mu(length(hat_mu));
% end

% for i=hat_Tau(length(hat_Tau)-1):hat_Tau(length(hat_Tau))
%     hat_s(i)=hat_mu(length(hat_mu));
% end;

end

