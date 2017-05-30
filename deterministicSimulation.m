function [Signal, Delta0, L0, Tau] = deterministicSimulation(N)


%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%  Firsty, we simulate a benchmark with 24 change points

    Tau=N*[0.1111 0.1213 0.1504 0.1707 0.2002 0.2203 0.2402 0.2511 0.2613 0.2712 0.3005 0.3302 0.3505 0.4012 0.5011 0.55003 0.5602 0.7004 0.7402 0.8004 0.8301 0.8504 0.9001 0.9501];
    delta_mu=[-1.1, 1, -1.25, 1.85, 1.7, 1, 2.2, -1.2, -2.5, 1.75, -1.5, -1.8, 1.1, 1.4, 1.5, -2.5, 1, -1, 1.25, -1.15, -2.2, 1.1, 1.4, -1.8 ]; 
  
    Delta0=min(abs(delta_mu));  % delta0 is the minimum jump size (in absolute value)
    Tau= round(Tau);   % Tau should be a vector of ingegers
        % Tau denotes the right configuration of change points on the mean

    K=length(Tau); % number of change points on the mean
    
%% Calculation of L0 the min distance between two successive change points
    
    for k=1:length(Tau)-1
        L(k)=Tau(k+1) - Tau(k);
    end
    L0=min(L);
    
 %% Simulation

    % we denote by Signal the right signal, that is the mean E(X)
    % Signal has jump of size delta_mu(k) at each change point Tau(k)
    %
    % on the first interval(i.e. before the first jump) the value of the
    % mean is mu0.

    mu0=4;  % value of means before the change time #1 
    
   clear Signal;
   Signal(1:N)=mu0;  %initialization
 
    mu=mu0;
    for k=1:length(Tau)-1
        mu=mu+delta_mu(k); % at each change point Tau(k) the mean jumps of delta_mu(k)
        for i=Tau(k)+1:Tau(k+1)
            Signal(i)=mu;  % 
        end;
    end;
    

end

