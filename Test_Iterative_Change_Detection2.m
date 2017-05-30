%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (May 2017)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data and function for change detection are in specific folders
close all;
clear all;

addpath('RessourcesFDpV');  % Folder with all the functions for change point detection

FigureNumber=1;   % number of the plotted figure
NEWSIMULATION=0;

%% 1) INPUT
if NEWSIMULATION==0
    load('/SimulationPaper.mat');
else
    %%  We simulate a benchmark
    SERIESLENGTH=30000;        % length of the data set
    SIGMA=2.7;        % standard deviation (std)

    %% Simulation of the observed signal
    [ObservedSeries, Signal, MinJump, MinLength, ChangePoints]=simulation2(SERIESLENGTH, SIGMA);  % with 24 change points
end;

%% Description of the variables X, Signal, Delta0, L0, Tau, K
% the function Simulation depends 
% on the lenght SERIESLENGTH of the series and the std SIGMA
% Output:
% ObservedSeries (X): a series of length SERIESLENGTH with std SIGMAa auround 
% Signal: the right signal, piecewise constant with change at times 
% ChangePoints (Tau): the position of the changes times
% MinLength (L0):  is the minimum distance between two successive change times 
% MinJump (Delta0): the minimum size of jump (in absolute value)

K=length(ChangePoints)-1; % number of change points on Signal
%  Warning : length(Tau)= K+1, since we add the last value Tau(K+1)=N


% we check the simulation, by plotting Figure 1 & 2.
FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid;
hold on;
set(gca,'LineWidth', 3);
plot(ObservedSeries,'y');
%plot(Signal,'r --','linewidth',2);
legend('\fontsize{20} Observed signal', 'Location','NE'); 
title(['\fontsize{24} The right number of  change points is  ', num2str(length(ChangePoints)-1)]);
hold off; 

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid;
hold on;
set(gca,'LineWidth', 3);
plot(ObservedSeries,'y');
plot(Signal,'r --','linewidth',2);
legend('\fontsize{20} Observed signal', '\fontsize{20} right signal', 'Location','NE'); 
title(['\fontsize{24} The right number of  change points is  ', num2str(length(ChangePoints)-1)]);
hold off; 


%% 2) we fix the extra-parameters of the iterative method
MinStudentConverge=52;     % we force a minimum distance between two successive change points 
                          % insuring convergence of Student law to normal law 
                          % and variance
               
KMAX=150;                   % maximum number of potential change points 
FALSE_ALARMRISK=0.025      % the risk for False Alarm
NON_DETECTIONRISK=0.01     % the risk level for Non Detection

TC0=5;          % the first threshold is 0.1*TC0, ie. 0.5
TCF=40;         % the maximum value of last threshold is 0.1*TCF, ie. 4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% STEP 0: Preliminary step 
%%
%% for estimating the std SIGMA assumed to be unknown
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0=60;   % window size for estimating variance SIGMA

[A, TresholdStep1, sigma2, UncertaintyTP, Dmin, MinStudentConverge, lambda2] = estimateParameters2(ObservedSeries, A0, MinJump, KMAX, FALSE_ALARMRISK, NON_DETECTIONRISK, MinStudentConverge);

% the function estimateParameters calculates the extra-parameters used for 
% determining the set of potential change points Tau1
% then cancelling the false alarm by iteratively increasing the threshold
% for t-values. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% STEP 1: The Filtered Derivative method
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1.1) Calculating FD function
% 
% Since we have changed the window size, we have again to calculate 
% Filtered Derivative function
%clear FD;
[FD]=filteredDerivative(ObservedSeries, A);

% Step 1.2: Selecting potential change points
%clear Tau0; 
%clear FD1;
[Step1PotentialChangePoints, FD1]=potentialChangePoints(FD, 1, SERIESLENGTH, Dmin, TresholdStep1, KMAX);

PotentialChangePoints=Step1PotentialChangePoints;  % Step1PotentialChangePoints 
                                                   % is the set of potential change points after Step 1 
 
                                                   
% % We calculate positive and negative threshold
% PositiveThresholds=TresholdStep1*ones(SERIESLENGTH);   % vector of positive threshold +TresholdStep1
% NegativeThresholds=-PositiveThresholds;         % vector of negative threshold -TresholdStep1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  STEP 2.2 : The Iterative procedure
%%
IndexCriticalTvalues=[TC0:TCF];

clear tval; clear hat_mu; clear hat_var; 
clear LengthPotentialSegments; clear SigD Dobs;
[Tvalues, MeanOnSegments,  VarianceOnSegments, LengthPotentialSegments, StdOnChangePoints, Dobs] = pvalues(ObservedSeries, PotentialChangePoints, UncertaintyTP);

%% OLD Version:
%clear Tau1s; 
% Tau1s=PotentialChangePoints;  % a temporary set of potential change points
%% It is useless to have a temporary vector of change points Tau1s

hatK(1)=length(PotentialChangePoints)-1;  % hatK is a vector or a real ???

%hatK=length(PotentialChangePoints)-1;

MinObservedTvalues= min(abs(Tvalues));  % Minimum Tvalue (in absolute value)
    
% tval1 and Tau1 should have the same length
Tvalues=[Tvalues NaN];

%% initialization

% M contains frames for the movie / animation
M(length(IndexCriticalTvalues)) = struct('cdata', [], 'colormap', []);
IndexMovie = 1; % index for Movie

% clear EstimateNumberTruePositiveTmp;
% clear MaxTvaluesFAtmp;
% clear hatFA; clear estim_FA;


if hatK>0
        MaxTvaluesFAtmp=abs(norminv((1-(1-FALSE_ALARMRISK)^(1/hatK))/2)); % Maximum Tvalues of False Alarms, 
end;                                                                   % Formula (3.10) in Hadouni-Dutheil-Bertrand (2017)

SetEstimateTruePositives=find(abs(Tvalues)>MaxTvaluesFAtmp);  % True Positive corresponds to t-values greater than MaxTvaluesFA
EstimateNumberTruePositiveTmp=length(SetEstimateTruePositives);  

EstimateNumberFATmp=hatK(1)-EstimateNumberTruePositiveTmp;        % Estimate number of false alarm

IndexCriticalTvalue=TC0;
EstimateNumberFA(IndexCriticalTvalue)=length(PotentialChangePoints);


while (min(Tvalues) < MaxTvaluesFAtmp && IndexCriticalTvalue<TCF)
    % we increment  IndexCriticalTvalue and calculate the critical t-value
    IndexCriticalTvalue=IndexCriticalTvalue+1;
    CriticalTvalue=0.1*IndexCriticalTvalue;

    % At step IndexCriticalTvalue, the set PotentialChangePoints is the
    % subset satisfying: abs(Tvalues(k)) >= CriticalTvalue
    [PotentialChangePoints ]= elimin(PotentialChangePoints, Tvalues, CriticalTvalue);
    % The function elimin, eliminates  PotentialChangePoints(k)
    %  for which  abs(Tvalues(k)) < CriticalTvalue
    [PotentialChangePoints]=[PotentialChangePoints SERIESLENGTH];
    PotentialChangePoints=round(PotentialChangePoints);
  
    % calculation of the t-values
    % clear Tvalues;
    [Tvalues, MeanOnSegments, Dobs] = pvalues(ObservedSeries, PotentialChangePoints, UncertaintyTP);
    
    % Tvalues and PotentialChangePoints should have the same length
    Tvalues=[Tvalues NaN];
    
    %% Estimation of Max tvalue False Alarm, Number of True Positives, and False Alarm
    % Since we are making a simulation, we know the True positives and the false alarm
    % in the set of potential change points.
    % But, we play the part of user of real dataset, who doesn't known the
    % rigth change points.
    
    EstimateNumberTruePositiveTmp=length(find(abs(Tvalues)>MaxTvaluesFAtmp));
    EstimateNumberFATmp=hatK(1)-EstimateNumberTruePositiveTmp;
    % we iterate one time the estimation of Max(tvalues False Alarms)
    if EstimateNumberFATmp>0
        MaxTvaluesFAtmp=abs(norminv((1-(1-FALSE_ALARMRISK)^(1/EstimateNumberFATmp))/2));
    else
        MaxTvaluesFAtmp=EstimateMaxTvaluesFA(IndexCriticalTvalue-1);
    end;
    % However,it is better to iterate two times
    EstimateNumberTruePositive(IndexCriticalTvalue)=length(find(abs(Tvalues)>MaxTvaluesFAtmp));
    EstimateNumberFA(IndexCriticalTvalue)=hatK(1)-EstimateNumberTruePositive(IndexCriticalTvalue);
    EstimateMaxTvaluesFA(IndexCriticalTvalue)=MaxTvaluesFAtmp;


    %% Separation of True Positives, False Alarms and computation of the corresponding t-values
    % Kowning the right configuration of changes, we can compute the number
    % of True positives and False alarm in the set of potential change
    % points, at each iteration.
    clear index_TruePositive; clear TrueEpsilon;
    [index_TruePositive  TrueEpsilon] = separationFA_TP(PotentialChangePoints,ChangePoints,A);
    
    NumberNonDetection(IndexCriticalTvalue)=length(ChangePoints)-sum(index_TruePositive);
    NumberFalseAlarm(IndexCriticalTvalue) =length(find(index_TruePositive<1)); %length(Tau0)-sum(index_TruePositive);
    NumberTruePositive(IndexCriticalTvalue)=sum(index_TruePositive)-1; % also =length(find(index_TruePositive>0))-1
    
    TvalueTPs=Tvalues(find(index_TruePositive>0));  %% set of t-value of true positives
    TvalueFAs=Tvalues(find(index_TruePositive<1));   % set of t-value of false alarms

    if NumberFalseAlarm(IndexCriticalTvalue)>0
        MaxTvalueFA(IndexCriticalTvalue)= max(TvalueFAs);
        theoretical_Max_tvalFA(IndexCriticalTvalue) = abs(norminv(FALSE_ALARMRISK/(2*NumberFalseAlarm(IndexCriticalTvalue))));
    end;
    MinTvalueTruePositives(IndexCriticalTvalue)= min(TvalueTPs);
    
%% Estimation of the shift at step it
%% Calculation N_k and Delta_k
    clear EstimateShifts; 
    EstimateShifts=deltak(PotentialChangePoints, UncertaintyTP, MinJump, sigma2);  %
    EstimateMinShift(IndexCriticalTvalue)=min(EstimateShifts);
    
    
%% Signal reconstruction and plotting
    clear EstimateSignal;  
    
    [EstimateSignal] = signalReconstruction(PotentialChangePoints, MeanOnSegments);
    
    %% computation of the average Square Error (ISE)   
    AverageSE(IndexCriticalTvalue)=sum((EstimateSignal-Signal).^2)/SERIESLENGTH;
    
    hatK(IndexCriticalTvalue)=length(PotentialChangePoints);
    MinObservedTvalues(IndexCriticalTvalue)= min(Tvalues);
    
    clear LengthSegments; 
    for k=2:length(PotentialChangePoints)
        LengthSegments(k)=PotentialChangePoints(k)-PotentialChangePoints(k-1);
    end;
    
    LengthSegments(1)=A;
    MinLengthSegments(IndexCriticalTvalue)=min(LengthSegments(2:length(LengthSegments)));
    
    
%% Plotting Signal recontruction
    figure(100+IndexCriticalTvalue)
    grid;
    hold on;
    plot(ObservedSeries,'y');
    plot(Signal,'g --','linewidth',2);
    plot(EstimateSignal, 'm', 'linewidth',2);
    legend('Observed signal X', 'right signal s',['estimated signal after Step ' num2str(i)], 'Location','NE'); 
    title(sprintf(['The true signal and the successive estimated  signals  for A=' num2str(A) ', C1=' num2str(TresholdStep1/MinJump) ' x delta_0',  '\n  with a ratio noise/signal= sigma/delta_0 = '  num2str(SIGMA/MinJump)]));
    u.Value = IndexMovie;
    M(IndexMovie)=getframe(gcf); 
    hold off;
    IndexMovie=IndexMovie+1;
 
    %% Plotting Histogram t-values TP and FA
    figure(200+IndexCriticalTvalue)
    grid;
    hold on;
    xlabel('t-value');
    h=histogram(TvalueTPs,20); 
    h.FaceColor = 'b'; %[0 0 1];  
    h.EdgeColor = 'r'; %[1 0 0];
    g=histogram(TvalueFAs);  %,20);
    g.FaceColor = [1 128/255 0];
    for k=0:7
        plot(MaxTvaluesFAtmp,k ,'+ r', 'linewidth',2);
    end;
    legend('Histogram of True Positive t-values', 'Histogram of False Alarms t-values', 'Location','NE'); 
    title(sprintf(['t-values  for A=' num2str(A) ', tc=' num2str(IndexCriticalTvalue/10)....
        ', C1=' num2str(TresholdStep1/MinJump) ' x delta_0',  '\n  with a ratio noise/signal SNR= sigma/delta_0 = '  num2str(SIGMA/MinJump)...
        ' , and L_0=' num2str(MinLength)]));
    hold off;
 %% Make a movie with the histogram TP and FA indexed by IndexCriticalTvalue 
    

end;

%% Viewing the movies
%pause;
%movie(gcf, M);

% movie(gcf, MovieHistograms);
%% Make a movie with the histogram (to be done)
    
NormL2TrueEpsilon=sqrt(mean(TrueEpsilon.^2))    

NormLinfinityTrueEpsilon=max(abs(TrueEpsilon))

Estimated_Uncertainty= UncertaintyTP

% %    norm(hatShift(it)-tcf/10)
% pd=makedist('Normal',0,1);
%  NonDetectionRisk= 1-cdf(pd,EstimateMinShift(IndexCriticalTvalue)-TCF/10) %^hat_NTP

TrueShift=MinJump/sigma2*sqrt((MinLength-2*UncertaintyTP)/2)
LastEstimatedShift=EstimateMinShift(IndexCriticalTvalue)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 3) Plotting the observed parameters
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid
ylabel('estimate numbers of change');
xlabel('iterative threshold (critical t-value)');
set(gca,'LineWidth', 3);
hold on;
plot(EstimateNumberTruePositive, 'g','linewidth',2);
plot(EstimateNumberFA, 'r','linewidth',2);
plot(NumberFalseAlarm, 'r --','linewidth',1);
plot(NumberNonDetection, 'm','linewidth',2);

legend('Estimate number of true positive','Estimate number of false alarm', 'residual number of false alarm', 'number of non detection');
%title(sprintf(['Histogram of the jumps size, mean jump=' num2str(mean_jump) ', std duration= '  num2str(std_jump) ' \n  number of change point= ' num2str(hat_Ki)]));  
hold off


FigureNumber=FigureNumber+1;
figure(FigureNumber);
    grid
    hold on;
    xlabel('iterative threshold (critical t-value)');
    ylabel('t-value');
 %   set(gca,'XTickLabel',idx_tc_set*0.1');
    plot(MinTvalueTruePositives, 'g --','linewidth',2);
    plot(MaxTvalueFA, 'r -.','linewidth',2);
    plot(EstimateMaxTvaluesFA, 'm -.','linewidth',2);
    plot(MinObservedTvalues, 'b --','linewidth',2);
    plot(NumberNonDetection, 'm','linewidth',2);
    plot(EstimateMinShift,'g -','linewidth',3);    
    legend('Min t-values for True Positives', 'Max t-values for False Alarms', ['Estimated Max tval FA at risk ' num2str(FALSE_ALARMRISK)],'minimum observed t-values','Number of Non Detection','Estimation of the Shift');  
    title(sprintf(['Min and max of t-values for True Positives and False Alarms \n  for A=' num2str(A)...
        ', L_0=' num2str(MinLength) ', C1=' num2str(TresholdStep1/MinJump)...
        ' x delta_0',  '\n  with a ratio noise/signal= sigma/delta_0 = '  num2str(SIGMA/MinJump)]));
    hold off;

FigureNumber=FigureNumber+1;
figure(FigureNumber);
hold on;
xlabel('local shift');
grid;
histogram(EstimateShifts,10);
title(sprintf(['Histogram of shifts at the end of iterative FDtV \n for A=' num2str(A) ', L_0=' num2str(MinLength) ', epsilon_0=' num2str(UncertaintyTP)]));
hold off;

% FigureNumber=FigureNumber+1;
% figure(FigureNumber);
% hold on;
% xlabel('Minimum Length of Segments');
% grid;
% plot(MinLengthSegments,'g --', 'linewidth',2);
% title(sprintf(['Minimum Length of Segments1 for iterative FDtV for A=' num2str(A) ...
%     '\n , where the right minimum length is L_0=' num2str(MinLength) ', epsilon_0=' num2str(UncertaintyTP)]));
% hold off;


FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid
hold on;
xlabel('\fontsize{36} iterative threshold (critical t-value) * 10');
ylabel('\fontsize{36} t-value');
set(gca,'LineWidth', 3);
%set(gca,'XTickLabel',idx_tc_set*0.1');
plot(MinTvalueTruePositives, 'g --','linewidth',3);
plot(MaxTvalueFA, 'r -.','linewidth',3);
plot(EstimateMaxTvaluesFA, 'm -.','linewidth',3);
plot(MinObservedTvalues, 'b --','linewidth',3);
plot(NumberNonDetection, 'm','linewidth',3);
legend('\fontsize{20} Min t-values for True Positives', '\fontsize{20} Max t-values for False Alarms', ['\fontsize{20} Estimated Max tval FA at risk ' num2str(FALSE_ALARMRISK)],'\fontsize{20} minimum observed t-values','\fontsize{20} Number of Non Detection');  
title(sprintf(['Min and max of t-values for True Positives and False Alarms \n  for A=' num2str(A)...
    ', L_0=' num2str(MinLength) ', C1=' num2str(TresholdStep1/MinJump) ' x delta_0',  '\n  with a ratio noise/signal= sigma/delta_0 = '  num2str(SIGMA/MinJump)]));
%     title('\fontsize{24} Min and max of t-values for True Positives and False Alarms for A=290, L_0=297, C1=0.37 with SNR=delta0/sigma=0.33');
hold off;



%% 3) COMPARISON with other methods
%
% 3.1) We compare with a non iterative FDpV procedure
% Recall:
% Step1PotentialChangePoints is the set of potential change points after Step 1

CriticalTvalue_FDpV=3.5;  

[Tvalues_FDpV, MeansOnSegments_FDpV, Dobs] = pvalues(ObservedSeries, Step1PotentialChangePoints, UncertaintyTP);
% Tvalues_FDp and PotentialChangePoints should have the same length
Tvalues_FDpV=[Tvalues_FDpV NaN];


%% signal reconstruction with FDpV non iterative 
[PotentialChangePoints]= elimin(Step1PotentialChangePoints, Tvalues_FDpV, CriticalTvalue_FDpV);
if(PotentialChangePoints(length(PotentialChangePoints)) ~= SERIESLENGTH)
    [PotentialChangePoints]=[PotentialChangePoints SERIESLENGTH];
end
PotentialChangePoints=round(PotentialChangePoints);

%clear tval; clear hat_mu; clear Dobs;
[Tvalues_FDpV, MeansOnSegments_FDpV, Dobs] = pvalues(ObservedSeries, PotentialChangePoints, UncertaintyTP);
% Tvalues_FDpV=[Tvalues_FDpV NaN];

%clear EstimatedSignal_FDpV;
[EstimatedSignal_FDpV] = signalReconstruction(PotentialChangePoints, MeansOnSegments_FDpV);

IntegralSquareError_iterative=AverageSE(IndexCriticalTvalue)
IntegralSquareError_FDpV=sum((EstimatedSignal_FDpV-Signal).^2)/SERIESLENGTH

 

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid;
hold on;
xlabel('\fontsize{36} \tau')
ylabel('\fontsize{36} \mu')
plot(ObservedSeries,'b');
plot(Signal,'r', 'linewidth',5);
plot(EstimateSignal, 'g', 'linewidth',4);
plot(EstimatedSignal_FDpV, 'm--', 'linewidth',3);
% plot(ObservedSeries,'y');
% plot(Signal,'g --','linewidth',2);
% plot(EstimateSignal, 'm', 'linewidth',2);
% plot(EstimatedSignal_FDpV, 'b --', 'linewidth',2);
legend('\fontsize{24} Observed signal X', '\fontsize{24} right signal s','\fontsize{24} estimated signal by iterative FDtV', '\fontsize{24} estimated signal by FDpV (non iterative)', 'Location','NE'); 
title(sprintf(['The right signal and the successive estimated  signals  for A=' num2str(A) ...
    ', L_0=' num2str(MinLength) ', C1=' num2str(TresholdStep1/MinJump) ' x delta_0', ...
    '\n  with an estimated signal-to-noise ratio SNR=delta_0/sigma=' num2str(MinJump/sigma2)]));
hold off;
    
    
%% 3.2)  We compare with FD method, with C2 chosen for having less than 0.01
%       False Alarms (Meisser et al.)

 load('q01.mat');
 load('Aset.mat');

TresholdFD=sigma2*q01(Aset((A-Aset(1))/10));

[PotentialChangePoints, FD1]=potentialChangePoints(FD, 1, SERIESLENGTH, Dmin, TresholdFD, KMAX);
[Tvalues, MeanOnSegments,  Dobs] = pvalues(ObservedSeries, PotentialChangePoints, UncertaintyTP);

[EstimatedSignalFD] = signalReconstruction(PotentialChangePoints, MeanOnSegments);

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid;
hold on;
xlabel('\fontsize{36} \tau')
ylabel('\fontsize{36} \mu')
plot(ObservedSeries,'b');
plot(Signal,'g', 'linewidth',5);
plot(EstimateSignal, 'r--', 'linewidth',4);
plot(EstimatedSignal_FDpV, 'm', 'linewidth',3);
plot(EstimatedSignalFD, 'c -.', 'linewidth', 2);
legend('\fontsize{20} Observed signal X', '\fontsize{20} right signal','\fontsize{20} estimated signal by iterative FDtV', '\fontsize{20} estimated signal by FDpV (non iterative)', '\fontsize{20} estimated signal by FD', 'Location','NE'); 
title(sprintf(['The true signal and the successive estimated  signals  for A=' num2str(A) ...
    ', L_0=' num2str(MinLength) ', C1=' num2str(TresholdStep1/MinJump) ' x delta_0', ...
    '\n  with an estimated ratio noise/signal= sigma/delta_0 = '  num2str(sigma2/MinJump)]));
hold off;
 
IntegralSquareErrorFD=sum((EstimatedSignalFD-Signal).^2)/SERIESLENGTH










% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
% %% Addendum A.1) Stat for Fred Dutheil
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [hat_Ki,jump, duration]=Stat_for_Fred(EstimateSignal,PotentialChangePoints);
% % 
% % mean_duration=mean(duration);
% % std_duration=std(duration);
% % 
% % Nb_rupture=hat_Ki
% % 
% % 
% % FigureNumber=FigureNumber+1;
% % figure(FigureNumber);
% % grid
% % xlabel('duration');
% % ylabel('Histogram');
% % set(gca,'LineWidth', 3);
% % hold on;
% % h=histogram(duration,30); 
% %     h.FaceColor = 'b'; %[0 0 1]; %Int?rieur des rectangles
% %     h.EdgeColor = 'r'; %[1 0 0];
% % title(sprintf(['Histogram of the time duration, mean duration=' num2str(mean_duration) ', std duration= '  num2str(std_duration) ' \n  number of change point= ' num2str(hat_Ki)]));  
% % hold off
% % 
% % 
% % % calculation delta_k 
% % abs_jump=abs(jump);
% % mean_jump=mean(abs(jump));
% % std_jump=std(abs(jump));
% % 
% % 
% % FigureNumber=FigureNumber+1;
% % figure(FigureNumber);
% % grid
% % xlabel('jump on Heart Rate (beats/minute)');
% % ylabel('Histogram');
% % set(gca,'LineWidth', 3);
% % hold on;
% % g=histogram(abs_jump,30); 
% %     g.FaceColor = 'b'; %[0 0 1]; %Int?rieur des rectangles
% %     g.EdgeColor = 'r'; %[1 0 0];
% % %legend('Instantaneous heart-rate (beat/minute)','piecewise constant mean heart-rate');
% % title(sprintf(['Histogram of the jumps size, mean jump=' num2str(mean_jump) ', std duration= '  num2str(std_jump) ' \n number of change point= ' num2str(hat_Ki)]));  
% % hold off
% %
% % FigureNumber=FigureNumber+1;
% % figure(FigureNumber);
% % grid
% % ylabel('jump on Heart Rate (beats/minute)');
% % xlabel('duration');
% % set(gca,'LineWidth', 3);
% % hold on;
% % plot( duration,jump, 'b +', 'linewidth',3);
% % %legend('Instantaneous heart-rate (beat/minute)','piecewise constant mean heart-rate');
% % title(sprintf(['Histogram of the jumps size, mean jump=' num2str(mean_jump) ', std duration= '  num2str(std_jump) ' \n  number of change point= ' num2str(hat_Ki)]));  
% % hold off
% % 
% %% END of STAT for Fred Dutheil
