function [windowIterative, lambda, nu]=choiceParametersIFDtV(SNR, SeriesLength, NonDetectionRisk, FalseAlarmRisk, MinStudentConverge, TC0)
%
%% We apply Proposition 3.9  in Bertrand & Hadouni(2017)
%
%% INPUT
% SeriesLength=30000;         % length of the series (dataset)
% NonDetectionRisk=0.02;    % risk of Non Detection 
% FalseAlarmRisk=0.05;      % risk of False Alarm
% MinStudentConverge=52;      % Minimum sample size insuring convergence of Student law to normal law 
%                             % and empirical variance to the theoretical one.       
% TC0=1;                      % first threshold for iterative FDtV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Monte-Carlo simulations provide the values nu for different SNRs

SNRs=[0.3  0.33  0.4  0.5 0.66 0.75  1  1.25 1.5 1.75  2   2.3  2.5   2.75  3   3.3];
Nus= [6.2  6.8    7   7.5  7.6  7.7  8   8   8.5  8.6  8.7  8.8  8.5  8.4  8.3  8.2];  %obtained by Monte-Carlo simulations

% We deduce nu for the observed SNR
Q=polyfit(SNRs,Nus, 3);
nu=polyval(Q,SNR);

%%%%%%%%%%%%%%%
%% 
%% Iterative FDtV
%%
MinWindowStudentConverge= MinStudentConverge +2*nu/(SNR*SNR);       
        % A2: minimum window size to insure convergence of Student law to Normal law

Kmax= floor(SeriesLength/MinWindowStudentConverge-1);           
        % the maximum number of potential change point satisfy: Kmax* A2 < N

xTmp=TC0+ norminv((1-NonDetectionRisk/2)^(2/Kmax));
mu=(norminv((1-NonDetectionRisk/2)^(1/Kmax)))/sqrt(xTmp*xTmp+nu);
lambda=1-mu;
 
MinWindowStep1=2*norminv((1-NonDetectionRisk/2)^(1/Kmax))^2/(mu*mu*SNR*SNR); 
        % A1: minimum window size for Step 1
MinWindowIterativeStep2=2*(xTmp*xTmp+nu)/(SNR*SNR);             
        % A3: minimum window size for Step 2

%%  we ajust Kmax= N/max(A1, A2,A3)
Kmax=floor(SeriesLength/max(max(MinWindowStep1,MinWindowStudentConverge), MinWindowIterativeStep2)-1);   
        % we ajuste Kmax by using   Kmax* Max(A1,A2,A3) < N


xTmp=TC0+ norminv((1-NonDetectionRisk/2)^(2/Kmax));  % an intermediate variable
                                                                           % used in definition of Mu, then Lambda
mu=(norminv((1-NonDetectionRisk/2)^(1/Kmax)))/sqrt(xTmp*xTmp+nu); 
                                    % Mu=1-Lambda, 
                                    % where
                                    % meaning of Lambda is given in the next lines
                                    
lambda=1-mu;    % the threshold of Step 1 is 
                                    % C1= Lambda * SNR * Sigma 
                                    % (Proposition 3.7 Hadouni-Dutheil-Bertrand, 2017)
MinWindowStep1=2*norminv((1-NonDetectionRisk/2)^(1/Kmax))^2/(mu*mu*SNR*SNR); 
                                                                           % minimum window size for Step 1
                                                                           % insuring zero non detection at RISK_NON_DETECTION/2
MinWindowIterativeStep2=2*(xTmp*xTmp+nu)/(SNR*SNR);   
                                    % minimum window size for Step 2
                                    % insuring zero non detection at RISK_NON_DETECTION/2 
                                    % and zero false alarm at RISK_FALSE_ALARM

windowIterative=max(max(MinWindowStep1,MinWindowStudentConverge), MinWindowIterativeStep2);                                    
end

