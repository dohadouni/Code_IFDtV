function [ IRS1 ] = incrementRatio(X)
% compute the ratio statistic of order 1
% corresponding to the
% INPUT
% X: a time series
% IRS: the associate Increment Ratio time series
 
N = length(X);

for k=1:(N-1)
    delta1(k)=X(k+1)-X(k);
end;

for k=1:(N-2)
    IRS1(k)=abs(delta1(k+1)+delta1(k))/(abs(delta1(k+1))+ abs(delta1(k)));
    if (delta1(k+1)==0 && delta1(k)==0)
        IRS1(k)=1;
    end;
end;



% epsilon=0.00001;
% 
% for k=1:(N-2)                  % Compute Increment Ratio time series
% IRS1(k)=abs(delta1(k+1)+delta1(k))/(epsilon+abs(delta1(k+1))+ abs(delta1(k)));
% %if (delta1(k+1)==0 && delta1(k)==0)
% %IRS1(k)=1;
% %end;
% end;

end

