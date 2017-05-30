function [EstimatedSignal, ChangePoints]=iterativeFDpV2(ObservedSeries, A0, MinJump, Kmax, NonDetectionRisk, FalseAlarmRisk)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%%  Pierre R. BERTRAND - Doha HADOUNI - Guillaume PAUGAM (April-September 2016)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the signal obtained by iterative FDpV method
%
% INPUT
% ObservedSeries (X): the observed series
% MinJump: the minimum jump size
% Kmax: maximum number of change points
% NonDetectionRisk: the risk level for non detection 
% FalseAlarmRisk: the risk level for false alarm
%
% OUTPUT
% EstimatedSignal: the signal after change detection 
% ChangePoints: the family of change times

%% The extra-parameters of the iterative method


minStudentConverge=32;  % we force a minimum distance between two successive change points 
                        % after substracting  intertitude on localization
                        % i.e. min(Tau0(k+1)-Tau0(k) - 2 * epsilon > MinStudentConverge
                        % insuring convergence of Student law to normal law and variance                 
tc0=10;     % the first threshold is 0.1*tc0=1
tcf=40;     % the maximum value of last threshold is 0.1*tcf=4


[A, thresholdStep1, sigma, uncertainty, Dmin, minStudentConverge] = estimateParameters3(ObservedSeries, A0, MinJump, Kmax, NonDetectionRisk, FalseAlarmRisk, minStudentConverge, tc0/10);
thresholdStep1
% % Previous version with estimateParameters2
% [A, thresholdStep1, sigma, uncertainty, Dmin, minStudentConverge] = estimateParameters3(ObservedSeries, A0, MinJump, Kmax, NonDetectionRisk, FalseAlarmRisk, minStudentConverge, tc0/10);

A=10*(floor(max(A,A0)/10)+1)

% we force the window size A to be greater than A0.
% Warning: in Test_Iterative_Change_Detection.m
% we need a window A being a multiple of 10

SNR=MinJump/sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% STEP 1: The Filtered Derivative method
%%
%%

% Step 1.1) Calculating FD function

clear FD;
[FD]=filteredDerivative(ObservedSeries, A);

% Step 1.2: Selecting potential change points (Step1= FilteredDerivative)
Dmin=A;
[changePointsStep1, FD1]=potentialChangePoints(FD, 1, length(ObservedSeries), Dmin, thresholdStep1, Kmax);
hatK=length(changePointsStep1)-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  STEP 2.2 : The Iterative procedure
%%

[Tvalues, MeanOnSegments, Dobs] = pvalues(ObservedSeries, changePointsStep1, uncertainty);
% tval1 and Tau1 should have the same length
Tvalues=[Tvalues NaN];

TvaluesStep1=Tvalues;

ChangePoints=changePointsStep1;
    
if length(ChangePoints)>1
        MaxTvalueFA=abs(norminv((1-(1-FalseAlarmRisk)^(1/hatK))/2)); 
end;
%% ELSE ?
%% what happen if length(changePointsStep1)=1?
%
setTruePositives=find(abs(Tvalues)>MaxTvalueFA);
NumberFalseAlarm=hatK-length(setTruePositives);

indexCriticalTvalue=tc0;

while (min(Tvalues) < MaxTvalueFA && indexCriticalTvalue<tcf) %%(min(tval) < hat_mtvFA+0.05 && it<tcf)
    indexCriticalTvalue=indexCriticalTvalue+1;
    criticalTvalue=0.1*indexCriticalTvalue;

    if max(Tvalues)>criticalTvalue
        [ChangePoints]= elimin(ChangePoints, Tvalues, criticalTvalue);
        [ChangePoints]=[ChangePoints length(ObservedSeries)];
        ChangePoints=round(ChangePoints);
%     else
%         changePoints=changePoints;
    end;

% calculation of the t-values
    [Tvalues, MeanOnSegments, Dobs] = pvalues(ObservedSeries, ChangePoints, uncertainty);
    % tval and Tau should have the same length
    Tvalues=[Tvalues NaN];
    
%% Estimation of Max tvalue False Alarm
    setTruePositives=find(abs(Tvalues)>MaxTvalueFA);
    NumberTruePositive=length(setTruePositives);
    NumberFalseAlarm=length(changePointsStep1)-NumberTruePositive-1;  %
    % we iterate the estimation of Max(tvalues False Alarms)
    if NumberFalseAlarm>0
        MaxTvalueFA=abs(norminv((1-(1-FalseAlarmRisk)^(1/NumberFalseAlarm))/2));       
%     else
%         MaxTvalueFA=estim_mtvFA(indexCriticalTvalue-1);
    end;
%    estim_mtvFA(indexCriticalTvalue)=MaxTvalueFA;
end;

[FalseAlarms, indexFA]=setdiff(changePointsStep1,ChangePoints);


%% Plotting T-values
figure(23)
grid;
hold on;
xlabel('\fontsize{36} \tau');
ylabel('\fontsize{36} t-value');
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'LineWidth', 3);
plot(changePointsStep1, TvaluesStep1, 'b *'); 
plot(ChangePoints, Tvalues, ' g o', 'linewidth',3);
plot(FalseAlarms, TvaluesStep1(indexFA), 'r o', 'linewidth',3); 
CriticalTvalues=MaxTvalueFA * ones(length(changePointsStep1));
plot(changePointsStep1,CriticalTvalues,'m --','linewidth',2);
legend('\fontsize{24} t-values of potential change points', '\fontsize{24} t-values of selected change points', '\fontsize{24} t-values of false alarms',['\fontsize{24} Threshold False alarm vs True Positive = ' num2str(MaxTvalueFA)],'location', 'north')
title('\fontsize{24} T-values of the change points for the signal estimated');
hold off;
    
% Signal reconstruction 
[EstimatedSignal] = signalReconstruction(ChangePoints, MeanOnSegments);
    % WARNING:
    % we should have length(MeanOnSegments)= length(changePoints)
end