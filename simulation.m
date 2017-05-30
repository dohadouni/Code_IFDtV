function [X, Signal, Delta0, L0, Tau]=simulation(N, sigma)
%
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




%%  Tau with 6 change points for the paper Hadouni-Dutheil-Bertrand (Figures 1-6)


Tau=N*[0.2004 0.2607 0.3506 0.4302 0.5012 0.7602];  % 6 change points
Tau= round(Tau);

K=length(Tau); % number of change points on the mean

% configuration of means
mu0=3;
delta_mu=[1.8, 1.85, 2, -1, 2.5, -1.75]; %values used in the paper Hadouna-Dutheil (2017)
    
    
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