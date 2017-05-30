%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (June 2016)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Monte-Carlo simulation for the extra parameters
%% We compute the number of false detection, the number of non detection and the unvertainty of localization
%% for different values of the extra-parameters A and SNR,
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% Data and function for change detection are in specific folders

addpath('RessourcesFDpV');  % Folder with all the functions for change point detection


%% We fix the extra-parameter

NUMBERMONTECARLOSIMULATION=250;   % Number of Monte-Carlo simulations
KMAX=150;           % maximum number of potential change points 
SERIESLENGTH=20000; % length of the series

SIGMASET=[1, 0.5, 2, 2.5];  % set of standard deviation (std)
WINDOWSET= [20:10:320];     % set of values used for window size A

FigureNumber=15;

%% 1) INPUT
%%  Firsty, we simulate a deterministic signal with 24 change points

[Signal, MinJump, MinLength, ChangePoints] = deterministicSimulation(SERIESLENGTH);
% DeterministicSimulation simulates a series: Signal of length SERIESLENGTH
% with change points at times ChangePoints (Tau), 
% mimimum jump size MinJump (Delta_0), 
% and mimimum length of segment MinLength (L0)

K=length(ChangePoints); % number of change points on the mean    
    
%% 2) Monte-Carlo simulation:
%% for each j=1:M, we simulate a sample,
%% then estimate the number of non detection and the number of false alarm
%% for different choices of the extra parameters A, and C1.

for Sigma=SIGMASET
    %% Calculation of Signal/Noise Ratio (dB)
    SNR=MinJump/Sigma;           % SNR is the Signal/Noise Ratio (where Signal=minimum jump size)
    SNRdB=log(SNR)/log(10)      % Signal/Noise Ratio (dB)
    
    for j=1:NUMBERMONTECARLOSIMULATION;
        %% simulating the noisy signal X = s+ rdn
        %  we add a Gaussian centered random variable with std Sigma  
        X=Signal+Sigma.*randn(1,SERIESLENGTH);
        for A=WINDOWSET
            %% Step 1: The Filtered Derivative map FD is calculated by a specific function:
            clear FD;
            [FD]=filteredDerivative(X, A);
            
            ThresholdStep1=min(abs(FD(ChangePoints)))*0.95 ;  % Threshold for Step1
                                                     % insuring zero (0) non detection
                                                     % of True Positives
                                                     
            Ratio1(A, j)=ThresholdStep1/MinJump;    % corresponding ratio
            
            %% Step 2: Selecting the set of potential change points (Tau0)
            clear PotentialChangePoints; clear FD1;
            [PotentialChangePoints, FD1]=potentialChangePoints(FD, 1, SERIESLENGTH, A, ThresholdStep1, KMAX);
            %
            % WARNING: the last value of Tau0 is N, which is not really a change point
%           %  clear HatK;
%           %  HatK=length(Tau0)-1;

            % we calculate the length of each segment
            clear hat_L;
            PotentialSegmentLengths(1)=PotentialChangePoints(1)-1;
            for k=1:length(PotentialChangePoints)-1;
                PotentialSegmentLengths(k+1)=PotentialChangePoints(k+1)-PotentialChangePoints(k);
            end;
            MinPotentialLength(A,j)=min(PotentialSegmentLengths);

            %%  Separation of the false alarm (FA) and right detection 
            %%
            %% We have Tau0 a vector of potential change points, 
            %% we want to extract two sub-vectors of False Alarms (FA) and True Positive (TP)
            %%
            clear index_TruePositive; 
            clear TruePositive;
            clear  NumberTP; clear TrueEpsilon;

            [index_TruePositive, TrueEpsilon] = separationFA_TP(PotentialChangePoints,ChangePoints,A);
            % we compute the indices of True Positives in the vector PotentialChangePoints
            % and the TrueEpsilon's = the error of position of TruePositive wrt ChangePoints

            NumberTruePositive=sum(index_TruePositive);  % we compute the number of true positives

            NumberNonDetection(A,j)=length(ChangePoints)-sum(index_TruePositive);  
                % Non detection are right change point which are not true positives
            NumberFalseAlarm(A,j)=length(PotentialChangePoints)-1-sum(index_TruePositive);
                % False alarm are potential change point which are not true positives
                % Warning: the last value Tau_{K+1} is in Tau0. Thus, we have to
                % substract one

           % Max_TrueEpsilon(A,j)=max(abs(TrueEpsilon));
          %  L1_TrueEpsilon(A)=mean(abs(TrueEpsilon));
            L2NormTrueEpsilon(A,j)=sqrt(mean(abs(TrueEpsilon).^2));

            lambda2(A, j)=0;
            NTP=0;
            if NumberTruePositive>0    
                for k=1:length(PotentialChangePoints)
                    if index_TruePositive(k)>0
                        NTP=NTP+1;
                        TruePositive(NTP)=PotentialChangePoints(k);
                    end;
                end;
                lambda2(A, j)= min(abs(FD(TruePositive)))/MinJump;
            end;
        end;
    end;

    %% Calculation of the mean value of the Monte Carlo calculations

    for A=WINDOWSET
        Mean_Ratio1(A) = mean(Ratio1(A, 1:NUMBERMONTECARLOSIMULATION));
        % std_Ratio1(A) = std(Ratio1(A, 1:M));
        % M2S_Ratio1(A) = Mean_Ratio1(A) -2.33* std(Ratio1(A, 1:M));  %std_Ratio1(A);

        Mean_lambda2(A) = mean(lambda2(A, 1:NUMBERMONTECARLOSIMULATION));
        % std_lambda2(A) = std(lambda2(A, 1:M));
        % at risk alpha=1%
        % M2S_lambda2(A) = Mean_lambda2(A) -2.33* std(lambda2(A, 1:M)); %std_lambda2(A);  

        MeanNumberNonDetection(A)= mean(NumberNonDetection(A,1:NUMBERMONTECARLOSIMULATION));
        MeanNumberFalseAlarm(A)= mean(NumberFalseAlarm(A,1:NUMBERMONTECARLOSIMULATION));
        % Std_FA(A)= std(NumberFalseAlarm(A,1:M));
        % M2S_FA(A)=MFA(A)+ 2.33*std(NumberFalseAlarm(A,1:M));  %Std_FA(A);

        Mean_hatL0(A)=mean(MinPotentialLength(A,1:NUMBERMONTECARLOSIMULATION));

        % Mean_LinfinyEpsilon(A)= mean(Max_TrueEpsilon(A,1:NUMBERMONTECARLOSIMULATION));
        % std_LinfinyEpsilon(A)= std(Max_TrueEpsilon(A,1:M));
        % M2S_LinfinyEpsilon(A)=  Mean_LinfinyEpsilon(A)+2.33*std(Max_TrueEpsilon(A,1:M));  %std_LinfinyEpsilon(A);

        MeanL2NormTrueEpsilon(A)= mean(L2NormTrueEpsilon(A,1:NUMBERMONTECARLOSIMULATION));
        stdL2NormTrueEpsilon(A)= std(L2NormTrueEpsilon(A,1:NUMBERMONTECARLOSIMULATION)); %std_L2Epsilon(A);
        M2SL2NormTrueEpsilon(A)= MeanL2NormTrueEpsilon(A)+2* std(L2NormTrueEpsilon(A,1:NUMBERMONTECARLOSIMULATION)); 
        % +2.33* stdL2NormTrueEpsilon
    end;


%     %% Statistical interpretation
% 
%     %% Aset_FD The set of admissible condition for Filtered Derivative method (without Step2)
% 
% 
%     %alpha=0.01;
%     %Aset_NDFA=find(MND(Aset)<0.05 & MFA(Aset)<0.05);
%     Aset_NDFA=find(MeanNumberNonDetection(WINDOWSET)<0.01 & MeanNumberFalseAlarm(WINDOWSET)<0.05);
%     FD_admissible_set=WINDOWSET(Aset_NDFA)
% 
%     if length(Aset_NDFA)>0
%         FD_admissible_set=WINDOWSET(Aset_NDFA);  %% The set of admissible values for FD method (without step 2)
%         Mean_change_point_square_error=MeanL2NormTrueEpsilon(WINDOWSET(Aset_NDFA))
%     %    At_risk_01=M2S_L2Epsilon(Aset(Aset_NDFA))
%     end;
% 
%     min(MeanNumberNonDetection(WINDOWSET))
% 
%     if min(MeanNumberNonDetection(WINDOWSET))<0.01
%         Aset_Step2=find(MeanNumberNonDetection(WINDOWSET)<0.01 & MeanNumberFalseAlarm(WINDOWSET)<50);
%     else
%         Aset_Step2=find(MeanNumberNonDetection(WINDOWSET)<0.05);
%         alpha=0.05
%     end;
% 
%     Step2_admissible_set=WINDOWSET(Aset_Step2)
% 
%     % if length(Aset_Step2)>0
%     %     M2S_L2_Epsilon=M2S_L2Epsilon(Aset(Aset_Step2))
%     % 
%     %     M2S_Linfinity_Epsilon=M2S_LinfinyEpsilon(Aset(Aset_Step2))
%     % end;
% 
% 
%     clear Mean_N0; clear Shift;
%      for i=WINDOWSET
%         Mean_N0(i)= Mean_hatL0(i) - 2*MeanL2NormTrueEpsilon(A); % M2S_L2Epsilon(i);   %Mean_LinfinyEpsilon(i);
%         Shift(i)=MinJump/Sigma*sqrt(Mean_N0(i)/2);
%     end;
% 
%     if length(Aset_Step2)>0
%         Mean_N0_STEP2=Mean_N0(WINDOWSET(Aset_Step2))
% 
%         Shift_H1_Step2=Shift(WINDOWSET(Aset_Step2))
% 
%     %    lambda2_at_risk01=M2S_lambda2(Aset(Aset_Step2))
%         MFA_Step2=MeanNumberFalseAlarm(WINDOWSET(Aset_Step2))
% 
%         True_Epsilon_divided_by_square_noise=MeanL2NormTrueEpsilon(WINDOWSET(Aset_Step2))/(Sigma/MinJump)^2
%         at_risk_01= M2SL2NormTrueEpsilon(WINDOWSET(Aset_Step2))/(Sigma/MinJump)^2
% 
% %        Max_True_Epsilon_divided_by_square_noise=Mean_LinfinyEpsilon(WINDOWSET(Aset_Step2))/(Sigma/MinJump)^2
% 
%         corresponding_Aset=WINDOWSET(Aset_Step2)
%     end;
% 
%     % real_lambda2=min(M2S_lambda2(Aset(Aset_Step2)))*0.95


    %% Plotting
    FigureNumber=FigureNumber+1;
    figure(FigureNumber)
    grid;
    hold on;
    xlabel('\fontsize{36} Window size A');
    set(gca, 'FontSize', 20, 'fontName','Times');
    axis([20  290  -2  65]);
    plot(WINDOWSET, MeanNumberFalseAlarm(WINDOWSET), 'r-.', 'linewidth',3);
    plot(WINDOWSET, MeanNumberNonDetection(WINDOWSET)*10, 'g','linewidth',3);
    plot(WINDOWSET, MeanL2NormTrueEpsilon(WINDOWSET), 'm -.', 'linewidth',3);
    clear y;
    y=0: 1 :  max(MeanNumberFalseAlarm(WINDOWSET));
    plot(MinLength, y, 'r +','linewidth',1);
    legend('\fontsize{24} Mean Number of False Alarm', '\fontsize{24} Mean Number of Non Detection *10', '\fontsize{24} Mean L2 norm of True Epsilon','location', 'NE')
    title(['\fontsize{30} Varying window size A with SNR =' num2str(SNRdB) 'dB (SNR= delta_0/sigma =' num2str(SNR) ') and L_0 = 198']);
    hold off;

end;