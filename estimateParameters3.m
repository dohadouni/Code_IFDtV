 
function [A, tresholdStep1, Sigma, Uncertainty, Dmin, N0, lambda] = estimateParameters3(ObservedSeries, A, MinJump, Kmax, NonDetectionRisk, FalseAlarmRisk, MinStudentConverge, tc0)


if nargin < 8
    tc1 = 5;
else 
    tc1  = tc0;
end
    

%% STEP 1: Estimate variance of signal
%
%  The variance of X over estimated the variance sigma^2
%  since it contains the jump and the rigth variance

Sigma=std(ObservedSeries);  % Temporary standard deviation

%% Calculation of Filtered Derivative function (FD)
%% and setection of potential change points

[FD]=filteredDerivative(ObservedSeries, A);  % The Filtered Derivative function is calculated by a specific function

Uncertainty=floor(10*(Sigma/MinJump)^2)+1;
% Since sigma over estimates sigma, epsilon over estimates the error of
% position
Dmin=2*Uncertainty + MinStudentConverge;
% we force a minimum distance between two successive change point

lambda=0.45;
tresholdStep1=lambda* MinJump; %The threshold is the minimum jump size (delta0) multiplided by a ratio lambda

clear Tau0;
[Tau0, FD1]=potentialChangePoints(FD, 1, length(ObservedSeries), Dmin, tresholdStep1, Kmax);
% the last value of Tau0 is N, which is not really a change point


%% Next we estimate the variance

[Tvalues MeanOnSegments  VarianceOnSegments LengthSegments StdOnChangePoints Dobs] = pvalues(ObservedSeries, Tau0, Uncertainty);

% We want to compute the estimated std, assumed to be constant.
% Old wrong formula
% sigmaTmp=sqrt(mean(VarianceOnSegments(1:length(VarianceOnSegments)-1)))

%  We need to cancell the last value
Sigma= sqrt(sum(LengthSegments(1:length(VarianceOnSegments)-1).*VarianceOnSegments(1:length(VarianceOnSegments)-1))/sum(LengthSegments(1:length(VarianceOnSegments))))

%% END OF VARIANCE ESTIMATION



%% 1) The following parameters can be automatically deduced from Sigma and MinJump.

SNR=MinJump/Sigma  % Signal/Noise Ratio


%% We apply Proposition 3.9 in Hadouni-Dutheil-Bertrand (2017)

[A, lambda, nu]=choiceParametersIFDtV(SNR, length(ObservedSeries), NonDetectionRisk, FalseAlarmRisk, MinStudentConverge, tc0)

tresholdStep1=lambda* MinJump;
Uncertainty=floor(nu/SNR^2)+1;
Dmin=A;
N0=MinStudentConverge;
end


