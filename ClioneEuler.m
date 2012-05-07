% fwd Euler for Clione swimming
clear all; hold off;

Total_Neurons = 8;  %Solve for this number of interacting Neurons
DT = 0.02;  %Time increment as fraction of time constant
Final_Time = 50;   %Final time value for calculation
Last = Final_Time/DT + 1;  %Last time step
Time = DT*[0:Last-1];  %Time vector
Tau = 0.8;  %Neural time constants in msec
TauR = 1.9;
TauSyn = 1.0;
ESyn = -0.92;   % for inhibitory synapses; 0 at excitatory
synThres=-0.20; %Threshold for IPSP conductance change
PIRoff=0;

k=9; % inhibitory gain
dorsalStim=0.5; % stimulus to dorsal neuron
I=dorsalStim; 

dorsalV=zeros(1,Last);
ventralV=zeros(1,Last);
dorsalV(1)=-0.70;
ventralV(1)=-0.70;
dorsalR=zeros(1,Last);
ventralR=zeros(1,Last);
dorsalR(1)=0.088;
ventralR(1)=0.088;
dorsalF=zeros(1,Last);
ventralF=zeros(1,Last);
dorsalF(1)=0;
ventralF(1)=0;
dorsalg=zeros(1,Last);
ventralg=zeros(1,Last);
dorsalg(1)=0;
ventralg(1)=0;

for t = 2:Last
	dorsalV(t)=dorsalV(t-1) + DT/Tau*(I - (17.81+47.71*dorsalV(t-1)+32.63*(dorsalV(t-1))^2)*(dorsalV(t-1)-0.55) - 26.0*dorsalR(t-1)*(dorsalV(t-1)+0.92) - k*(dorsalg(t-1))*(dorsalV(t-1)-ESyn)); 
    ventralV(t)=ventralV(t-1) + DT/Tau*(0 - (17.81+47.71*ventralV(t-1)+32.63*(ventralV(t-1))^2)*(ventralV(t-1)-0.55) - 26.0*ventralR(t-1)*(ventralV(t-1)+0.92) - k*(ventralg(t-1))*(ventralV(t-1)-ESyn)); 
    if t==150
        I=0;
    end
    
    dorsalR(t)=dorsalR(t-1) + DT/TauR*(1.03+1.35*(dorsalV(t-1))-dorsalR(t-1)+(dorsalV(t-1)+0.38)^2*3.3*PIRoff);
    ventralR(t)=ventralR(t-1) + DT/TauR*(1.03+1.35*(ventralV(t-1))-ventralR(t-1)+(ventralV(t-1)+0.38)^2*3.3*PIRoff);
    
    if ventralV(t-1)-synThres > 0
        Hstep=1;
    else
        Hstep=0;
    end
    dorsalF(t)=dorsalF(t-1) + DT/TauSyn*(-dorsalF(t-1) + Hstep);
    if dorsalV(t-1)-synThres > 0
        Hstep=1;
    else
        Hstep=0;
    end
    ventralF(t)=ventralF(t-1) + DT/TauSyn*(-ventralF(t-1) + Hstep);
    
    dorsalg(t)=dorsalg(t-1) + DT/TauSyn*(-dorsalg(t-1) + dorsalF(t-1));
    ventralg(t)= ventralg(t-1) + DT/TauSyn*(-ventralg(t-1) + ventralF(t-1));
end;


figure(1), ZA = plot(Time, 100*dorsalV, 'r', Time, 100*ventralV-150, 'b-'); set(ZA, 'LineWidth', 2);
ylabel('V (mV'); xlabel('Time (ms)');

VV = -0.9:0.01:1.5;
DVdt = -0.5*((1.37 + 3.67*VV + 2.51*VV.^2).*(VV - 0.55) - dorsalStim/13)./(VV + 0.92);
DRdt = 1.35*VV + 1.03;
figure(2), ZB = plot(VV, DVdt, 'k-', VV, DRdt, 'b-', dorsalV, dorsalR, 'r-'); axis([-1, 0.6, 0, 1]);
set(ZB, 'LineWidth', 2); axis square;
% %Next lines calculate spike rate
% Spikes = (X(1, 1:Last - 1) < -0.2).*(X(1, 2:Last) >= -0.2);
% SpkTime = zeros(1, sum(Spikes));
% Nspk = 1;  %Number of spike
% for T = 1:length(Spikes);  %Calculate spike rate for all interspike intervals
% 	if Spikes(T) == 1; SpkTime(Nspk) = T*DT; Nspk = Nspk + 1; end;
% end;
% Final = length(SpkTime);
% Rates = 1000./(SpkTime(2:Final) - SpkTime(1:Final - 1));
% Leng = length(Rates);
% Red_Rate = mean(Rates(Leng/2:Leng))
% %Next lines calculate blue spike rate
% Spikes = (X(3, 1:Last - 1) < -0.2).*(X(3, 2:Last) >= -0.2);
% SpkTime = zeros(1, sum(Spikes));
% Nspk = 1;  %Number of spike
% for T = 1:length(Spikes);  %Calculate spike rate for all interspike intervals
% 	if Spikes(T) == 1; SpkTime(Nspk) = T*DT; Nspk = Nspk + 1; end;
% end;
% Final = length(SpkTime);
% BRates = 1000./(SpkTime(2:Final) - SpkTime(1:Final - 1));
% Leng = length(BRates);
% Blue_Rate = mean(BRates(Leng/2:Leng))
