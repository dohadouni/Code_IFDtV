function [index_TruePositive  TrueEpsilon] = separationFA_TP(hatTau,Tau,A)
% DESCRIPTION
% separate the false alarms and the true positive 
% in the vector Tau0 of potential change points
% Tau0(j) is a true positive iff it is at a distance smaller than A from a
% true change point Tau(k)

% INPUT: 
% A:      window size
% Tau:    the true configuration of change times
% hatTau: the estimate configuration of change times

% OUPUT
% index_TruePositive: position of True Positive in hatTau
% TrueEpsilon: values of the error between each true positive and the corresponding
% change point


index_TruePositive(1:length(hatTau))=0;  %we initialize at 0 the presence of a true positive
                              % at rank j in hatTau 
TrueEpsilon(1:length(hatTau))= NaN;                          


    
% Next, we test wether a potential change point is in the vicinity of each
% change point Tau(k). 
% Vicinity = in the interval(Tau(k)-A, Tau(k)+A)
for k=1:length(Tau)
    M=find((hatTau>Tau(k)-A));
    % M1=find((hatTau(M)<Tau(k)+A));
    if ( hatTau(M(1))< Tau(k)+A )  
        index=M(1);

        %    M(2) can also be at a distance of Tau(k) smaller than A
        % if ( hatTau(M(2))< Tau(k)+A )  
        if length(M)>1
            if (abs(hatTau(M(2))-Tau(k))<abs(hatTau(M(1))-Tau(k)))
                index=M(2);
            end;
        end;

        TrueEpsilon(index)= hatTau(index)-Tau(k);

        index_TruePositive(index)=1;  % the potential change point # k is a true positive
    end; 
end;
setEpsilon=find(index_TruePositive>0);
TrueEpsilon=TrueEpsilon(setEpsilon);
end



