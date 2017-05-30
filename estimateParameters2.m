 
function [A, tresholdStep1, Sigma, Uncertainty, Dmin, N0, lambda] = estimateParameters2(ObservedSeries, A, MinJump, Kmax, alpha1, alpha2, N0, tc0)



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
Dmin=2*Uncertainty +N0;
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

%% 1) The following parameters can be automatically deduced from sigma/delta0
noises = [0.70   0.90  1.00    1.20  1.50  1.70 2.00  2.30  2.50  2.80   3.00];
lambda_obs =[0.56  0.56  0.524   0.53  0.53  0.53  0.48  0.45  0.44  0.41   0.37];
nu_obs =  [7.90  7.7   7.6  7.4  7.2  7.3 6.8  6.1  5.6  4.9  4.6];

P=polyfit(noises,lambda_obs,2);
Q=polyfit(noises,nu_obs,2);


%% 2) Choice of the extra-parameters from the estimated variance

SNR=MinJump/Sigma;  % Signal/Noise Ratio

lambda=polyval(P,1/SNR);     % the threshold for Filtered Derivative (Step1) is C1=lambda2*delta0 
nu=polyval(Q,1/SNR);        % the localization error is estimated by nu*(noise/signal)^2=nu*(sigma/delta0)^2

tresholdStep1=lambda* MinJump;


%% NEW VERSION
% modification calculation NO

hatK=length(Tau0)-1;
if hatK>0
        hat_mtvFA=abs(norminv((1-(1-alpha2)^(1/hatK))/2)); 
end;
setTP=find(abs(Tvalues)>hat_mtvFA);
hat_NTP=length(setTP);

init_shift = tc1 / 10;
Shift = init_shift + norminv((1-alpha1));

if hat_NTP>0
    Shift = init_shift + norminv((1-alpha1)^(1/hat_NTP));
end;
%Shift=1+ norminv(1-alpha1); 
N0Tmp=N0;
N0=2*Shift^2/SNR^2; % Modification of a wrong formula
N0=max(N0Tmp, floor(N0)+1);

% New version
Uncertainty=floor(nu/SNR^2)+1;

%% The window size can be automatically determined

A=10*floor((N0+2*Uncertainty)/10);
% Warning: in Test_Iterative_Change_Detection.m
% we need a window A being a multiple of 10

Dmin=A;

end

