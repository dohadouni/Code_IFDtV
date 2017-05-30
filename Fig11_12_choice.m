%% Choice of the extra-parameters 
%%
%% for FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (May 2017)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% We apply Formulas in Bertrand & Hadouni(2017)
%% Propositions 3.7 and Proposition 3.9
%%

clear all;
close all;

addpath('RessourcesFDpV');  % Folder with all the functions for change point detection

FigureNumber=100;   % number of the plotted figure

SERIESLENGTH=30000;         % length of the series (dataset)
RISK_NON_DETECTION=0.02;    % risk of Non Detection 
RISK_FALSE_ALARM=0.05;      % risk of False Alarm
MINSTUDENTCONVERGE=32;      % Minimum sample size insuring convergence of Student law to normal law 
                            % and empirical variance to the theoretical one.       
TC0=1;                      % first threshold for iterative FDtV

% Extra-parameters are calculate as a function of Signal/Noise Ratio (SNR)

SNRs=[0.3  0.33  0.4  0.5 0.66 0.75  1  1.25 1.5 1.75  2   2.3  2.5   2.75  3   3.3];
NUs= [6.2  6.8    7   7.5  7.6  7.7  8   8   8.5  8.6  8.7  8.8  8.5  8.4  8.3  8.2];  %obtained by Monte-Carlo simulations
%
% which can be translate in decibels
SNRdBs=log(SNRs)/log(10);

SNR=1/2.9701
SNRdB=log(SNR)/log(10)
[AIterative, lambda, nu]=choiceParametersIFDtV(SNR, SERIESLENGTH, RISK_NON_DETECTION, RISK_FALSE_ALARM, MINSTUDENTCONVERGE, TC0)

A2= MINSTUDENTCONVERGE +2.*NUs./(SNRs.*SNRs);       % minimum window size to insure convergence of Student law to Normal law

Kmax= floor(SERIESLENGTH./A2-1);           % the maximum number of potential change point satisfy: Kmax* A2 < N
Xtemp=abs(norminv((1-(1-RISK_FALSE_ALARM).^(2./Kmax))/2))+ norminv((1-RISK_NON_DETECTION/2).^(2./Kmax)); % an temporary vector
                                                                     % used for the definition of Mu and Lambda

Mu=(norminv((1-RISK_NON_DETECTION/2).^(1./Kmax)))./sqrt(Xtemp.*Xtemp+NUs);  % Mu=1-Lambda, 
                                                                     % where
                                                                     % the meaning of Lambda is given in the next lines

Lambda=1-Mu;                                                         % the threshold of Step 1 is 
                                                                     % C1= Lambda * SNR * Sigma 
                                                                     % (Proposition 3.7 Hadouni-Dutheil-Bertrand, 2017)

A1=2*norminv((1-RISK_NON_DETECTION/2).^(1./Kmax)).^2./(Mu.*Mu.*SNRs.*SNRs);  % minimum window size for Step 1
                                                                           % insuring zero non detection at RISK_NON_DETECTION/2
                                                            
A3=2*(Xtemp.*Xtemp+NUs)./(SNRs.*SNRs);                                              % minimum window size for Step 2
                                                                           % insuring zero non detection at RISK_NON_DETECTION/2 
                                                                           % and zero false alarm at RISK_FALSE_ALARM

Tc_FDpV=abs(norminv((1-(1-RISK_FALSE_ALARM).^(2./Kmax))/2));               % critical threshold of t-values for FDpV


%%%%%%%%%%%%%%%%%%%%%%%%
%% Enhanced bound

%% 1) we ajust Kmax=N/ max(A1, A2,A3)
Kmax_ajusted=floor(SERIESLENGTH./max(max(A1,A2), A3));    % actually, the maximum number of potential change point satisfy: 
                                        % Kmax* Max(A1,A2,A3) < N

% Next we computed Xx, Mu, Lambda, A1, A3, Tc 
% with the ajusted values  Kmax1.
Xx1=abs(norminv((1-(1-RISK_FALSE_ALARM).^(2./Kmax_ajusted))/2))+ norminv((1-RISK_NON_DETECTION/2).^(2./Kmax_ajusted));
Mu1=(norminv((1-RISK_NON_DETECTION/2).^(1./Kmax_ajusted)))./sqrt(Xx1.*Xx1+NUs);
Lambda1=1-Mu1;

A1_ajusted=2*norminv((1-RISK_NON_DETECTION/2).^(1./Kmax_ajusted)).^2./(Mu.*Mu.*SNRs.*SNRs);  % minimum window size for Step 1
A3_ajusted=2*(Xx1.*Xx1+NUs)./(SNRs.*SNRs);                                             % minimum window size for Step 2

Tc_ajusted=abs(norminv((1-(1-RISK_FALSE_ALARM).^(2./Kmax_ajusted))/2));                    % critical threshold of t-values for FDpV

%% End of ajusted Kmax and deduced variables


%%%%%%%%%%%%%%%
%% 
%% Iterative FDtV
%%

Xx_iterative=TC0+ norminv((1-RISK_NON_DETECTION/2).^(2./Kmax));
Mu_iterative=(norminv((1-RISK_NON_DETECTION/2).^(1./Kmax)))./sqrt(Xx_iterative.*Xx_iterative+NUs);
Lambda_iterative=1-Mu_iterative;
 
A1_iterative=2*norminv((1-RISK_NON_DETECTION/2).^(1./Kmax)).^2./(Mu_iterative.*Mu_iterative.*SNRs.*SNRs); 
                                                                        % minimum window size for Step 1
A3_iterative=2*(Xx_iterative.*Xx_iterative+NUs)./(SNRs.*SNRs);             % minimum window size for Step 2

%%  we ajust Kmax= N/max(A1, A2,A3)
Kmax_iterative=floor(SERIESLENGTH./max(max(A1_iterative,A2), A3_iterative));   % we ajuste Kmax by using  
                                                                    % Kmax* Max(A1,A2,A3) < N


Xx_iterative=TC0+ norminv((1-RISK_NON_DETECTION/2).^(2./Kmax_iterative));  % an intermediate variable
                                                                           % used in definition of Mu, then Lambda
Mu_iterative=(norminv((1-RISK_NON_DETECTION/2).^(1./Kmax_iterative)))./sqrt(Xx_iterative.*Xx_iterative+NUs); 
                                    % Mu=1-Lambda, 
                                    % where
                                    % meaning of Lambda is given in the next lines
                                    
Lambda_iterative=1-Mu_iterative;    % the threshold of Step 1 is 
                                    % C1= Lambda * SNR * Sigma 
                                    % (Proposition 3.7 Hadouni-Dutheil-Bertrand, 2017)
A1_iterative=2*norminv((1-RISK_NON_DETECTION/2).^(1./Kmax_iterative)).^2./(Mu_iterative.*Mu_iterative.*SNRs.*SNRs); 
                                                                           % minimum window size for Step 1
                                                                           % insuring zero non detection at RISK_NON_DETECTION/2
A3_iterative=2*(Xx_iterative.*Xx_iterative+NUs)./(SNRs.*SNRs);   
                                    % minimum window size for Step 2
                                    % insuring zero non detection at RISK_NON_DETECTION/2 
                                    % and zero false alarm at RISK_FALSE_ALARM

%% Calculation by a function


for k=1:length(SNRs)
    SNR=SNRs(k);
[AIts(k), lambdaIt(k), nuIt(k)]=choiceParametersIFDtV(SNR, SERIESLENGTH, RISK_NON_DETECTION, RISK_FALSE_ALARM, MINSTUDENTCONVERGE, TC0);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Plotting
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ajusted Kmax allows decreases the minimum window size A1 and A3
% for zero nondetection at risk RISK_NON_DETECTION 
% and zero false alarm at risk RISK_FALSE_ALARM
% as shown by the following figure

FigureNumber=FigureNumber+1;
figure(FigureNumber);
hold on;
grid;
xlabel('\fontsize{20} log_{10}(SNR)=log_{10} (Signal/Noise) (dB)');
ylabel('\fontsize{20} Minimum window size A');
set(gca,'LineWidth', 3);
plot(SNRdBs, Kmax_ajusted, 'b','linewidth',2);
plot(SNRdBs, A1, 'b +','linewidth',2);
plot(SNRdBs, A2, 'g -.','linewidth',2);
plot(SNRdBs, A3, 'r','linewidth',2);
plot(SNRdBs, A3_ajusted, 'm','linewidth',2); plot(SNRdBs, A3_ajusted, 'b +','linewidth',2);
plot(SNRdBs, Kmax_ajusted, 'r +','linewidth',2);
legend('\fontsize{20} Kmax ajusted','\fontsize{20} minimum window size A1 for zero non detection at Step 1 (FD)', '\fontsize{20} minimum window size A2 for long enough segments (N0>=32)', '\fontsize{20} minimum window size A3 for exactness at step 2 (FDpV)','\fontsize{20} minimum window size A1=A3 for FDpV with Kmax ajusted');
%legend('\fontsize{24} minimum window size A1 for zero non detection at Step 1 (FD)', '\fontsize{24} minimum window size A3 for exactness at step 2 (FDpV)','\fontsize{24} minimum window size A2 for long enough segments (N0>32)','\fontsize{24} A min for iFDpV', '\fontsize{24} minimum window size A1_iterative for Step 1 of iFDtV');
title(sprintf(['Domain (value of the window in function of SNR) \n where FDpV works \n for n=', num2str(SERIESLENGTH) , ',  alpha_1=', num2str(RISK_NON_DETECTION) ,' and alpha_2=', num2str(RISK_FALSE_ALARM),'\n enhanced version on Kmax']));
hold off;


% 
% Fig_number=Fig_number+1;
% figure(Fig_number);
% hold on;
% grid;
% xlabel('log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
% ylabel('x, tc, lambda*10');
% set(gca,'LineWidth', 3);
% plot(SNRdB, Xx, 'b','linewidth',2);
% plot(SNRdB, 10*Lambda, 'g -.','linewidth',2);
% plot(SNRdB, Tc_FDpV, 'r','linewidth',2);
% % plot(SNR, A4, 'm','linewidth',2);
% legend('lower bound of the Shift for FDpV', 'lambda *10', 'Tc: critical t-value for FDpV'); 
% title(sprintf(['parameters of FDpV method where FDpV works \n for n=', num2str(N) , ',  alpha_1=', num2str(RISK_NON_DETECTION) ,' and alpha_2=', num2str(RISK_FALSE_ALARM)]));
% hold off;

%% Ajusted Kmax has almost no effect of the threshold for step 1
% C1= Lambda * SNR *Sigma, where Lambda is plotted in the next figure
% but decreases the threshold of Step 2 for FDpV Tc_FDpV, 
% and the shift between True Positives t-values and False alarm ones.
%% Note that Lambda is almost constant and equal to 0.50.

FigureNumber=FigureNumber+1;
figure(FigureNumber);
hold on;
grid;
xlabel('log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
ylabel('x, tc, lambda*10');
set(gca,'LineWidth', 3);
plot(SNRdBs, Xtemp, 'b','linewidth',2);
plot(SNRdBs, 10*Lambda, 'g ','linewidth',2);
plot(SNRdBs, Tc_FDpV,  'r -.','linewidth',2);

plot(SNRdBs, Xx1, 'b +','linewidth',3);
plot(SNRdBs, 10*Lambda1, 'r + ','linewidth',3);
plot(SNRdBs, Tc_ajusted,'m +','linewidth',3);
 
plot(SNRdBs, Xx1, 'b --','linewidth',1);
plot(SNRdBs, 10*Lambda1, 'r -- ','linewidth',1);
plot(SNRdBs, Tc_ajusted,'m --','linewidth',1);

legend('lower bound of the Shift for FDpV', 'lambda *10', 'Tc: critical t-value for FDpV','lower bound of the Shift for FDpV (Kmax enhanced)', 'lambda *10  (Kmax enhanced)', 'Tc: critical t-value for FDpV  (Kmax enhanced)'); 
title(sprintf(['Parameters of FDpV method where FDpV works \n for n=', num2str(SERIESLENGTH) , ',  alpha_1=', num2str(RISK_NON_DETECTION) ,' and alpha_2=', num2str(RISK_FALSE_ALARM)]));
hold off;




%% Figures in the paper Hadouni-Dutheil, Bertrand, 2017
%% with the same numbers

FigureNumber=10;
FigureNumber=FigureNumber+1;
figure(FigureNumber);
% figure(11); %same number than in the paper
hold on;
grid;
xlabel('\fontsize{30} log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
ylabel('\fontsize{30} A min');
set(gca,'LineWidth', 3);
plot(SNRdBs, A1_ajusted, 'b *','linewidth',3);
plot(SNRdBs, A3_ajusted, 'r','linewidth',3);
plot(SNRdBs, A2, 'g','linewidth',3);
plot(SNRdBs, A3_iterative, 'm +','linewidth',3); 
plot(SNRdBs, A1_iterative, 'b','linewidth',3);
legend('\fontsize{24} minimum window size A1 for zero non detection at Step 1 (FD)', '\fontsize{24} minimum window size A3 for exactness at step 2 (FDpV)','\fontsize{24} minimum window size A2 for long enough segments (N0>=32)','\fontsize{24} A min for iFDpV', '\fontsize{24} minimum window size A1 iterative for Step 1 of iFDtV');
title([' \fontsize{20} Domain (value of the window in function of SNR) where FDpV works for n=', num2str(SERIESLENGTH) ,',  alpha_1=', num2str(RISK_NON_DETECTION) ,', and alpha_2=', num2str(RISK_FALSE_ALARM)]);
hold off;

%FigureNumber=10;
FigureNumber=FigureNumber+1;
figure(FigureNumber);
% figure(11); %same number than in the paper
hold on;
grid;
xlabel('\fontsize{30} log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
ylabel('\fontsize{30} A min');
set(gca,'LineWidth', 3);
plot(SNRdBs, A1_ajusted, 'b *','linewidth',3);
plot(SNRdBs, A3_ajusted, 'r','linewidth',3);
plot(SNRdBs, A2, 'g','linewidth',3);
plot(SNRdBs, A3_iterative, 'm +','linewidth',3); 
plot(SNRdBs, A1_iterative, 'b','linewidth',3);
plot(SNRdBs, AIts, '-- r','linewidth',3);
plot(SNRdBs, AIts, 'o r','linewidth',3);
legend('\fontsize{20} minimum window size A1 for zero non detection at Step 1 (FD)', '\fontsize{20} minimum window size A3 for exactness at step 2 (FDpV)','\fontsize{20} minimum window size A2 for long enough segments','\fontsize{20} A min for iFDpV', '\fontsize{20} minimum window size A1 iterative for Step 1 of iFDtV','\fontsize{20} minimum window size for iFDtV');
title([' \fontsize{20} Domain (value of the window in function of SNR) where FDpV works for n=', num2str(SERIESLENGTH) ,',  alpha_1=', num2str(RISK_NON_DETECTION) ,', and alpha_2=', num2str(RISK_FALSE_ALARM)]);
hold off;

FigureNumber=FigureNumber+1;
figure(FigureNumber);
% figure(12); %same number than in the paper
hold on;
grid;
xlabel('\fontsize{30} log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
ylabel('\fontsize{30}  x, tc, lambda*10');
set(gca,'LineWidth', 3);
plot(SNRdBs, Xx1, 'b','linewidth',3);
plot(SNRdBs, 10*Lambda1, 'g ','linewidth',3);
plot(SNRdBs, Tc_ajusted,  'r ','linewidth',3);
% iFDtV
plot(SNRdBs, Xx_iterative, 'b +','linewidth',3);
plot(SNRdBs, 10*Lambda_iterative, 'g + ','linewidth',3);
plot(SNRdBs, Xx_iterative, 'b --','linewidth',2);
plot(SNRdBs, 10*Lambda_iterative, 'g --','linewidth',2);
legend('\fontsize{24} Lower bound of the Shift for FDpV', '\fontsize{24} lambda *10', '\fontsize{24} Tc: critical t-value for FDpV','\fontsize{24} Shift for iterative FDtV','\fontsize{24} lambda *10 for iterative FDtV' ); 
title(['\fontsize{20} Parameters of FDpV method where FDpV works for n=', num2str(SERIESLENGTH) , ', alpha_1=', num2str(RISK_NON_DETECTION) ,' and alpha_2=', num2str(RISK_FALSE_ALARM)]);
hold off;

 
