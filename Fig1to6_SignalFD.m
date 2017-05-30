%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Change point analysis
%%  by iterative FDpV method
%%
%% by Pierre R. BERTRAND (80%), Doha HADOUNI (10%) & Guillaume PAUGAM (10%)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Les variables globales constantes en full majuscules 
% Les variables globales qui peuvent varier (pas l ideal en matlab cest  pour l oriente object) avec une majuscule au debut
% les fonctions commencent par une minuscule (et NON une majuscule comme en C++) 
% les variables locales en minuscules

%% This program plots Filtered Derivative function without noise 
%% with K=6 change points, with minimum jump Delta0, 
%% and minimum distance between two successive change times L0 (or minimum length)
%% for different window size A<L0/2, L0/2<A<L0, L0<A and A=1.5 * L0 > L0.


%% 1) INPUT
%%  Firsty, we simulate a benchmark

close all;
clear all;

addpath('RessourcesFDpV');  % Folder with all the functions for change point detection

SERIESLENGTH=2000;  % length of the series



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% STEP 1: The Filtered Derivative method
%%
%%

%Step1.0 (Preliminaries)
%       We fix the extra-parameters of Filtered Derivative and 1st selection of potential change points


KMAX=50;        % maximum number of potential change points 
LAMBDA=0.45;    % ratio for calculating the threshold C1

% with or without noise:
Sigmas=[0, 1];

% number of the plotted figure
FigureNumber=0;   

for Sigma=Sigmas
    %% Simulation of the observed signal

    [ObservedSeries, Signal , MinJump, MinLength, ChangePoints]=simulation(SERIESLENGTH, Sigma);

    %% Description of the variables X, Signal, Delta0, L0, Tau, K
    % the function Simulation depends 
    % on the lenght N of the series and the std Sigma
    % Output:
    % X: a series of length N with std Sigma auround 
    % Signal: the right signal, piecewise constant with change at times 
    % Tau: the position of the 6 changes times
    % L0:  is the minimum distance between two successive change times or
    % minimum length
    % Delta0: the minimum size of jump (in absolute value)

    K=length(ChangePoints)-1; % number of change points on the mean
    %  Warning : length(Tau)= K+1, since we add the last value Tau(K+1)=N
    
    C1=LAMBDA*MinJump;   %The threshold is delta0 multiplided by a ratio LAMBDA

    % We calculate positive and negative threshold
    PositiveThresholds=C1*ones(SERIESLENGTH);   % vector of positive threshold +C1
    NegativeThresholds=-PositiveThresholds;     % vector of negative threshold -C1
    
    % window size for computation of Filtered Derivative function

    Windows=[round(MinLength/2), MinLength-10, MinLength+30, round(MinLength*1.5) ]
    % Windows(1)= round(L0/2);   % a window size such that 2*A<L0
    % Windows(2)=L0-10;          % a window size such that A<L0<2*A 
    % Windows(3)=L0+30;          % a window size such that A>L0
    % Windows(4)=round(L0*1.5);  % a window size such that A>>L0.
    
    
    for A=Windows
        Dmin=A;         %minimum distance between two successive potential change points. 
                        % Set equal to the window size A.

        % Step 1.1) Calculating the filtered derivative function (FD)
        %
        clear FD;
        [FD]=filteredDerivative(ObservedSeries, A);

        % Step 1.2: Selecting potential change points
        clear Tau0; 
        clear FD1;
        [Tau0, FD1]=potentialChangePoints(FD, 1, SERIESLENGTH, Dmin, C1, KMAX);


        % Step 1.3: Plotting

        FigureNumber=FigureNumber+1;
        figure(FigureNumber);
        grid;
        hold on;
        set(gca, 'FontSize', 20, 'fontName','Times');
        plot(ObservedSeries,'b','linewidth',4);
        plot(Signal,'r--', 'linewidth',4);
        plot(FD,'m','linewidth',3);
        plot(PositiveThresholds,'g','linewidth',3);
        plot(NegativeThresholds,'g','linewidth',3);
        for k=1:length(Tau0)-1
           clear y;
           y=0: .05 : Signal(Tau0(k));
            plot(Tau0(k), y,'. r ');
        end;
        xlabel('\fontsize{36} \tau')
        ylabel('\fontsize{36} \mu')
        legend('\fontsize{26} Observed signal', '\fontsize{26} Right signal', '\fontsize{26} FD function','\fontsize{26} threshold +/-C_1', 'Location','NE');
        title(sprintf(['The filtered derivative function (FD)  with extra-parameters: window size A=' num2str(A) ' and threshold C_1=' num2str(C1) '\n for  the signal with noise sigma=' num2str(Sigma) ', minimum jump delta_0=' num2str(MinJump) ', and minimum length  L_0 = ' num2str(MinLength) ]));
        hold off; 
    end;
end; 