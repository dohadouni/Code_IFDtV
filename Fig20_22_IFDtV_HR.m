%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%) (May 2017)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Function used
%
% iterativeFDpV(X, A, delta0, Kmax, alphaND, alphaFA);
% 


close all;
clear all;

FigureNumber=20;   % number of the plotted figure

%% Data and function for change detection are in specific folders

addpath('RessourcesFDpV');  % Folder with all the functions for change point detection

%% 1) INPUT
load('Paris06_R1.mat');
%  load('Y1_Marseille.mat');  % uncomment to have a series of 110,000
%  heartbeats
RRs=Xf./1000;  %Xf is the RR series in milisecond. We translate in second.
RRsTmp=RRs(1: length(RRs)-950);  % the last minutes of Marathon HR are aberrant
clear RRs;
RRs=RRsTmp;

% For heartbeat
MinHRJump=10   %the minimine size of change on HeartRate

% Next, we can calculate  the minimine size of jumps on RR (in second)
MinRRJump= mean(RRs)^2/60* MinHRJump;

KMAX=500;       % maximum number of potential change points 
FALSE_ALARMRISK=0.025    % the risk level for max tvalFA
NON_DETECTIONRISK=0.01    % the risk level for Non Detection

A0=180;
%%% ITERATIVE FDpV
tic
[PieceWiseRRs, ChangePoints]=iterativeFDpV2(RRs, A0, MinRRJump, KMAX, NON_DETECTIONRISK, FALSE_ALARMRISK);
toc


%% OUTPUT= picts
%% Traduction in times (in minutes) and Heart Rate (beat/mn)

% Traduction of RR in HR in beats/minute
HR=60./RRs;
HRs=60./PieceWiseRRs;

%% traduction in times
clear times;

times(1)=RRs(1); %times in second
for i=1:length(RRs)-1
    times(i+1)=times(i)+RRs(i+1);
end;
times=times./60;  %time in minutes

% Calculation time resolution 
for k=1:length(ChangePoints)-1
    Lengths(k)=ChangePoints(k+1)-ChangePoints(k);
    timeDurations(k)=times(ChangePoints(k+1))-times(ChangePoints(k));
end;
MinLength=min(Lengths)
MinDuration=min(timeDurations);

MeanDuration=mean(timeDurations);
StdDuration=std(timeDurations);
DurationHours= length(RRs)*mean(RRs)/3600
NumberChange=length(ChangePoints)-1

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('\fontsize{36} Time in minute');
ylabel('\fontsize{36} Heart Rate in beat/mn');
set(gca,'LineWidth', 3);
hold on;
plot(times,HR,'b');
title('\fontsize{36} Instantaneous heart-rate during 6.26 hours');
hold off

FigureNumber=FigureNumber+1;
figure(FigureNumber);
grid
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('\fontsize{36} Time in minute');
ylabel('\fontsize{36} Heart Rate in beat/mn');
set(gca,'LineWidth', 3);
hold on;
plot(times,HR,'b');
plot(times, HRs, 'r', 'linewidth',3);
legend('\fontsize{24} Instantaneous heart-rate (beat/minute)','\fontsize{24} Piecewise constant mean heart-rate');
title(['\fontsize{25} Instantaneous heart-rate during' num2str(DurationHours) 'hours and its piecewise constant mean']);
% for A=' num2str(A) ' and C_1='0.0124');
hold off


FigureNumber=FigureNumber+2;
figure(FigureNumber);
grid
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('\fontsize{36} Time in minute');
ylabel('\fontsize{36} Heart Rate in beat/mn');
set(gca,'LineWidth', 3);
hold on;
set(gca,'LineWidth', 3);
plot(times, HRs, 'r', 'linewidth',4);
title(['\fontsize{36} Piecewise constant mean heart-rate with ' num2str(NumberChange) ' change points']);
hold off

