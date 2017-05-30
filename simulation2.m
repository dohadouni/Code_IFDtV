function [X, Signal, Delta0, L0, Tau]=simulation(N, sigma)

% Simulate a signal s corresponding to the configuration of change point Tau, and the observation X  
%   
% INPUT
% N: number of data
% sigma: noise
%
% OUTPUT
%
% X: the observed signal, of length N and std sigma
% Signal: the unknown signal of length N and without noise
% Delta0: the minimum size of jumps on the mean
% L0: the minimum distance between two successive change times
% Tau: the configuration of change times.


%%  Tau with 24 change points

Tau=N*[0.1111 0.1213 0.1504 0.1707 0.2002 0.2203 0.2402 0.2511 0.2613 0.2712 0.3005 0.3302 0.3505 0.4012 0.5011 0.5503 0.5602 0.7004 0.7402 0.8004 0.8301 0.8504 0.9001 0.9501  1];
Tau= round(Tau);

K=length(Tau); % number of change points on the mean

% configuration of means
    mu0=4;
    
    delta_mu=[-1.1, 1, -1.25, 1.85, 1.7, 1, 2.2, -1.2, -2.5, 1.75, -1.5, -1.8, 1.1, 1.4, 1.5, -2.5, 1, -1, 1.25, -1.15, -2.2, 1.1, 1.4, -1.8 ]; 
   %delta0= 1; 
    
%% Calculation of L0 the min distance between two successive change points 
%% and delta0 minimum jump size (in absolute value)
    
Delta0=min(abs(delta_mu)); 
if length(Tau)>1
    for k=1:length(Tau)-1
        L(k)=Tau(k+1) - Tau(k);
    end
    L0=min(L)
else
    L0=[];
end


%% Simulation

% we denote by Signal the right signal, that is the mean E(X)
% Signal has jump of size delta_mu(k) at each change point Tau(k)

Signal(1:N)=mu0;


mu=mu0;
for k=1:length(Tau)-1
    mu=mu+delta_mu(k);
    Signal(Tau(k)+1:Tau(k+1))=mu;
end;     
   
% next we add a Gaussian centered random variable with std sigma    
    
Z=randn(1,N);
X=Signal+sigma.*Z;    
    

end