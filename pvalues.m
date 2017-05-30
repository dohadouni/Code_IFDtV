function [Tvalues MeanOnSegments  VarianceOnSegments LengthSegments StdOnChangePoints Dobs] = pvalues(ObservedSeries, PotentialChangePoints, Uncertainty)

%% INPUT
% ObservedSeries(X):   the observed time series
% PotentialChangePoints (Tau0): set of potential change points  
% Uncertainty (Epsilon0) : uncertainty on localisation of change points
% Warning: The condition 2*epsilon<A should be satistied
%
%% OUTPUT
% MeanOnSegments: empirical mean of X on each segment
% VarianceOnSegments: empirical variance of X on each segment
% LengthSegments: length of the segments
% StdOnChangePoints: sqrt(S'_1^2/L1 +S'_2^/L2) 
%          empirical stds ont two successive segments araund a potential change point
% Tvalues: absolute value of tvalue for all potential change points = Dobs/SigD
%
%% Local variables
%
% tau1: a temporary vector of change points
% k: index (integer)


Tvalues =[];
MeanOnSegments  =[];
VarianceOnSegments =[];
LengthSegments =[];
StdOnChangePoints =[];
Dobs=[];

%% We have the condition  Tau(k+1)-Tau(k)> 2* epsilon  for k=0:length(PotentialChangePoints)
%% with the convention Tau(0)=0 and Tau(length(PotentialChangePoints)+1) = length(ObservedSeries).
%% We build a vector Tau1 satisfying these conditions
%% with 2*epsilon<A and Tau(k+1)-Tau(k)> 2*A.

clear tau1;
tau1(1)=1;
for k=2:(length(PotentialChangePoints)+1)
        tau1(k)=PotentialChangePoints(k-1); 
end;
tau1(length(PotentialChangePoints)+2)=length(ObservedSeries);

%% We caclulate the mean on each boxes after cancelling a vicinity of size
%% 2*epsilon0
% clear MeanOnSegments;
 
for k=1:(length(tau1)-1)
    MeanOnSegments(k)=mean(ObservedSeries(tau1(k)+Uncertainty:tau1(k+1)-Uncertainty));
end

%% We compute Dobs(k) which is the difference between mean on the right and
%% the left

for k=1:length(PotentialChangePoints) 
    Dobs(k)=abs(MeanOnSegments(k+1)-MeanOnSegments(k));
end;

%% We compute Y, defined by Y(j)=X(j)-mu(j) and Y(j)=0 if abs(j-Tau(k))<epsilon
for k=2:(length(tau1)-1)
    Y(tau1(k)+Uncertainty:tau1(k+1)-Uncertainty)=ObservedSeries(tau1(k)+Uncertainty:tau1(k+1)-Uncertainty)-MeanOnSegments(k);
end

%% we have a case disjonction for k=1 and k=Ntau1
Y(tau1(1):tau1(1)+Uncertainty-1)=0;
Y(tau1(1)+Uncertainty:tau1(2)-Uncertainty)=ObservedSeries(tau1(1)+Uncertainty:tau1(2)-Uncertainty)-MeanOnSegments(1);

Y(tau1(length(tau1))-Uncertainty+1:tau1(length(tau1)))=0;
Y(tau1(length(tau1)-1)+Uncertainty:tau1(length(tau1))-Uncertainty)=ObservedSeries(tau1(length(tau1)-1)+Uncertainty:tau1(length(tau1))-Uncertainty)-MeanOnSegments(length(tau1)-1);


%% We compute the length of the segment and
%% the local std  for all potential change point

for k=1:length(tau1)-1
    LengthSegments(k)=tau1(k+1)-tau1(k)-2*Uncertainty+1;
    VarianceOnSegments(k) = var(Y(tau1(k)+Uncertainty:tau1(k+1)-Uncertainty));
end;

for k=1:length(PotentialChangePoints)-1
    StdOnChangePoints(k)=sqrt(VarianceOnSegments(k)/LengthSegments(k)+  VarianceOnSegments(k+1)/LengthSegments(k+1));
    Tvalues(k)=Dobs(k)/StdOnChangePoints(k);
end;

%%%% nargout addendum for backwards compatibility with old pvalues2.m file.
if nargout == 3
	VarianceOnSegments = Dobs;
end

end


