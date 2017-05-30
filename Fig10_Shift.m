%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (September 2016)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The number of false alarm and the number of true positive are smaller 
%% than the maximum number of potential change points


clear all;
close all;

KMAX=350;      % maximum number of potential change points 
NON_DETECTIONRISK=0.02;    % risk of Non Detection 
FALSE_ALARMRISK=0.05;      % risk of False Alarm


TruePositiveNumbers=[1:KMAX];
FalseAlarmNumbers=KMAX-TruePositiveNumbers;
MaxFalseAlarmTvalues=abs(norminv((1-(1-FALSE_ALARMRISK).^(1./FalseAlarmNumbers))/2));
MinTruePositiveTvalues=norminv((1-NON_DETECTIONRISK).^(1./TruePositiveNumbers));
Shifts=MaxFalseAlarmTvalues+MinTruePositiveTvalues;

S0=max(Shifts)

figure(10)   % number of the figure in the article
hold on;
grid;
xlabel('\fontsize{26} number of true positives');
ylabel('\fontsize{26} Shift');
set(gca,'LineWidth', 3);
plot(MaxFalseAlarmTvalues, 'b','linewidth',2);
plot(MinTruePositiveTvalues, 'g','linewidth',2);
plot(Shifts, 'r','linewidth',2);
legend('\fontsize{26}\Psi^{-1}','\fontsize{26} \Phi^{-1}', '\fontsize{26} Shift', 'Location','NE');
title(['\fontsize{26} The functions \Phi^{-1} and \Psi^{-1} and the Shift for Kmax=', num2str(KMAX) ]);
hold off;

