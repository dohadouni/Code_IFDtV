function [K,jump, duration]=Stat_for_Fred(s,Tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   STAT  for Fred
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation time resolution 

K=length(Tau)-1;

for k=1:length(Tau)-1
    duration(k)=Tau(k+1)-Tau(k);
end;


% calculation delta_k in beats/minute

for k=1:length(Tau)-1
    jump(k)=s(Tau(k)+2) -s(Tau(k)-2);
end;
end
