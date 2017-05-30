function [Tau1 ]= elimin(Tau0, Tvalues, Threshold_Tvalues)
%% Description of the function  
%  Eliminate the potential change point Tau0(k) 
%  for which the t-value abs(Tvalue(k)) <Threshold_Tvalues
Tau1=[];
numberTau1=0;
for k=1:length(Tvalues)
    if abs(Tvalues(k))> Threshold_Tvalues
        numberTau1=numberTau1+1;
        Tau1(numberTau1)=Tau0(k);
    end
end

end

