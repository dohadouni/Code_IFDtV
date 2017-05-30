function [FD]=filteredDerivative(X, A)
%%
%% Description of the function
%% compute A times the filtered derivative (FD) function deduced from the time series X
%% input X, with a window size A.
% INPUT:
% X: vector= time series X(1),... X(N);
% A: window size
%
% OUPUT
% DF: Filtered Derivative function (vector of length N)
%
% Local variables:
% N: length of the series
% AFD= FD *A

N=length(X);

% We note FD the Filtered Derivative function (vector or time series of length N)
% We note AFD the Filtered Derivative function multiplied by A (the window size)
% We can recursively computate the function AFD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization to vector 0 
 AFD = zeros(N,1);
 %
 % AFD(A)=sum(X(A+1:2*A))-sum(X(1:A)); % the first value corresponds to  k=A
 % Iterative computation of AFD(k) for A <= k < N-A
 for k=A:N-A-1
         AFD(k+1)=AFD(k)+X(k+A+1)-2*X(k+1)+X(k-A+1); % recursive computation of function AFD   
 end
% division by A 
FD=AFD/A;

end
