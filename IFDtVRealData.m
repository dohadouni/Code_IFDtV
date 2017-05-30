%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (May 2017)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%% Data and function for change detection are in specific folders
addpath('RessourcesFDpV');  % Folder with all the functions for change point detection

%% Function used
% iterativeFDpV(X, A, delta0, Kmax, alphaND, alphaFA);
KMAX=500;       % maximum number of potential change points 
FALSE_ALARMRISK=0.05    % the risk level for max tvalFA
NON_DETECTIONRISK=0.01    % the risk level for Non Detection

A0=60;

FigureNumber=0;   % number of the plotted figure

%% 1) INPUT
load('ObservedSeries.mat')


MinJump=1   %the minimine size of change on HeartRate


%%% ITERATIVE FDpV
tic
[PieceWiseMeans, ChangePoints]=iterativeFDpV2(ObservedSeries, A0, MinJump, KMAX, NON_DETECTIONRISK, FALSE_ALARMRISK);
toc


%% OUTPUT= Plotting

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('\fontsize{36} Time');
ylabel('\fontsize{36} Values of the series');
set(gca,'LineWidth', 3);
hold on;
plot(ObservedSeries,'b');
plot(PieceWiseMeans, 'r', 'linewidth',3);
legend('\fontsize{24} Instantaneous series','\fontsize{24} Piecewise constant mean');
title('\fontsize{25} Instantaneous series and its piecewise constant mean');
hold off
