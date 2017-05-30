%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%%  Shift between t-values of false alarm and true positives
%%  and impact on distribution of abslolute value of t-values for true positives and false alarms
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (May 2017)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;

SNRs=[0.3  0.33  0.4  0.5 0.66 0.75  1  1.25 1.5 1.75  2   2.3  2.5   2.75  3   3.3];    % Set of Signal/Noise Ratio
NUs= [6.2  6.8    7   7.5  7.6  7.7  8   8   8.5  8.6  8.7  8.8  8.5  8.4  8.3  8.2];    % Set of coefficient NU 
                                                                                         % obtained by Monte-Carlo simulations

SNRdBs=log(SNRs)/log(10);   %  Signal Noise Ratio in decibel (dB)

TruePositiveUncertainties=floor(NUs .* SNRs.^(-2))+1;    % Formula (3.4) Hadouni-Dutheil-Bertrand, 2017


figure(7)  % same numerotation as in the paper
hold on;
grid;
plot(SNRdBs, TruePositiveUncertainties, 'r', 'linewidth', 3);
plot(SNRdBs, NUs, 'g', 'linewidth', 3);
plot(SNRdBs, TruePositiveUncertainties, 'b +', 'linewidth', 3);
plot(SNRdBs, NUs, 'm +', 'linewidth', 3);
xlabel('\fontsize{24} log_{10}(SNR)=log_{10} Signal/ Noise (dB)');
ylabel('\fontsize{24} \epsilon_0(SNR) and \nu(SNR)');
set(gca,'LineWidth', 3);
legend('\fontsize{24} Epsilon(SNR)','\fontsize{24} \nu(SNR)');
title(['\fontsize{24} The parameters Epsilon and nu depending on Signal/Noise Ratio (SNR)' ]); %\n for alpha_1=', num2str(alphaND) ,' and alpha_2=', num2str(alphaFA)]));
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Distribution of t-value for False Alarms (FA) 
%% and True Positives (TP) with a shift Delta0=4 or 7.5
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAMPLESIZE=10000;        % number of simulations
BINNUMBER=30;   % Number of bins for the histogram

% number of the plotted figure
FigureNumber=7;

Shifts=[4, 6, 7.5];    % shift for t-values of True Positives

FalseAlarmTvalues=randn(1,SAMPLESIZE);     % T-values of False Alarms

for Delta0=Shifts;              
    TruePositiveTvalues=FalseAlarmTvalues+Delta0; % T-values of True Positives
    % Plotting
    FigureNumber=FigureNumber+1;
    figure(FigureNumber)
        grid;
        hold on;
        xlabel('\fontsize{26} |t-value|');
        ylabel('\fontsize{26} Number of Potential Change Points');
        h=histogram(abs(TruePositiveTvalues),BINNUMBER);    % We are plotting absolute value of t-values
        h.FaceColor = 'b'; %[0 0 1] True Positives are blue
        h.EdgeColor = 'r'; %[1 0 0] False Alarms are orange
        g=histogram(abs(FalseAlarmTvalues), BINNUMBER);      % for both False Alarms and True Positives
        g.FaceColor = [1 128/255 0];
        legend('\fontsize{24} Histogram of True Positive t-values', '\fontsize{24} Histogram of False Alarms t-values', 'Location','NE'); 
        title(['\fontsize{30} Distritution of abolute values of t-values for a shift Delta0=' num2str(Delta0)]);
        hold off;
end;    
