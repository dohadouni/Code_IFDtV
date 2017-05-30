function [Tau0, FD1]=potentialChangePoints(FD, T0, Tf, Dmin, ThresholdStep1, Kmax)
%% Description of the function
% we detected the potential change point, corresponding to a absolute
% value of DF larger than the threshold ThresholdStep1, with a maximum number of
% change point Kmax.
% we use abs(FD)>Threshold <==> abs(AFD) > A*Threshold

%% INPUT
% FD:  Filtered Derivative function with a window size A
% T0: first time
% Tf: last time (or final time)
% Dmin: minimum distance between two successive potential change points
    % for example Dmin=A: window size
% Threshold: threshold
% Kmax: maximum number of potential change point

% WARNING:
% FD is defined for times T0+A<= t <=Tf-A

%% OUTPUT
% Tau: vector of potential change point
% FD1: the new Filtered Derivative function 
%      ater having set to 0 the vicinities of size A of each potential change point
%%
clear N; clear absDF; clear Tau;

FD1=FD;

N=length(FD);
absFD=zeros(N,1);

for i=(T0+Dmin):(Tf-Dmin)
    absFD(i)=abs(FD(i)); % absFD= abs(FD)  except outside (t0, tf)
end;

%%%%%%%%%%%%%%%%%%%% Initialisation: Searching potential change point #1
k=1; 
[Dmax(k) Tau(k)]=max(absFD);

%% Searching following change points
while Dmax(k)>ThresholdStep1 && (k<=Kmax-2) 
       for i=max(Tau(k)-Dmin,1):(min(Tau(k)+Dmin,N))
           absFD(i)=0; % we cancel FD(i) around Tau(k) on a window of size +/-Dmin
           FD1(i)=0;
       end;     
   k=k+1;
   [Dmax(k) Tau(k)]=max(absFD);
end;
Tau(k)=Tf;  %The last value is overwritten 
% since the last time is a change point

% we sort the potential change points
Tau0=sort(Tau);

end