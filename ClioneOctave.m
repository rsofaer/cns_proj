% fwd Euler for Clione swimming

function XT=Clione (kin,dorsalStimin,doPlot)
  global DT Final_Time Last Time Tau TauR TauSyn ESyn synThres PIRoff I k;
  DT = 0.01;  %Time increment as fraction of time constant
  Final_Time = 50;   %Final time value for calculation
  Last = Final_Time/DT + 1;  %Last time step
  Time = DT*[0:Last-1];  %Time vector
  Tau = 0.8;  %Neural time constants in msec
  TauR = 1.9;
  TauSyn = 1.0;
  ESyn = -0.92;   % for inhibitory synapses; 0 at excitatory
  synThres=-0.20; %Threshold for IPSP conductance change
  PIRoff=0;

  %Set argument defaults
  if(nargin<3)
    doPlot = 1;
  end
  if(nargin>=2)
    dorsalStim = dorsalStimin;
  end
  if(nargin>=1)
    k = kin;
  end
  if(nargin <= 1)
    dorsalStim = 0.5; % stimulus to dorsal neuron
  end
  if(nargin == 0)
    k = 9; % inhibitory gain
  end

  I = [zeros(1,150).+dorsalStim zeros(1,Last-150)];
  I = [I; zeros(1,Last)];

  Xinit = [
    -0.70; %dorsalV
    -0.70; %ventralV
    0.088; %dorsalR
    0.088; %ventralR
    0; %dorsalF
    0; %ventralF
    0; %dorsalG
    0; %ventralG
  ];

  Xtitles = [ 
  "dorsalV"
  "ventralV"
  "dorsalR"
  "ventralR"
  "dorsalF"
  "ventralF"
  "dorsalG"
  "ventralG"];
  X = [Xinit zeros(8,Last-1)];


  % Forward Euler:
  for t = 2:Last
    X(:,t) = X(:,t-1) + derivative(X(:,t-1),t-1).*DT;
  end;

  %The next two lines use Octave's built-in solver.
  %[X,istate,msg] = lsode(@derivative,Xinit,Time);
  %X = X';

  if(doPlot)
    figure(1), ZA = plot(Time, X(1,:), 'r', Time, X(2,:), 'b-'); set(ZA, 'LineWidth', 2);
    ylabel('V (mV'); xlabel('Time (ms)');

    VV = -0.9:0.01:1.5;
    DVdt = -0.5*((1.37 + 3.67*VV + 2.51*VV.^2).*(VV - 0.55) - dorsalStim/13)./(VV + 0.92);
    DRdt = 1.35*VV + 1.03;
    figure(2), ZB = plot(VV, DVdt, 'k-', VV, DRdt, 'b-', X(1,:), X(3,:), 'r-'); axis([-1, 0.6, 0, 1]);
    set(ZB, 'LineWidth', 2); axis square;
  end
  XT = [Time; X];
end
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

function xprime = derivative(Xt,t)
  t = t+1;
  xprime = [
    Vderivative(t,Xt,1,3,7);
    Vderivative(t,Xt,2,4,8);
    Rderivative(t,Xt,3,1);
    Rderivative(t,Xt,4,2);
    Fderivative(t,Xt,5,2);
    Fderivative(t,Xt,6,1);
    Gderivative(t,Xt,7,5);
    Gderivative(t,Xt,8,6)
    ];
end

%Dorsal V' means you want dorsal R, dorsal G, dorsal V
function xvprime = Vderivative(t,Xt,V_index,R_index,G_index) %All indicies are for the same side.
  global I ESyn Tau k;
  xvprime=(1/Tau)*(I(V_index,floor(t)) - (17.81+47.71*Xt(V_index)+32.63*(Xt(V_index))^2)*(Xt(V_index)-0.55) - 26.0*Xt(R_index)*(Xt(V_index)+0.92) - k*(Xt(G_index))*(Xt(V_index)-ESyn)); 
end

%Dorsal R' means dorsal V
function xrprime = Rderivative(t,Xt,R_index, V_index)
  global TauR PIRoff;
  xrprime = (1/TauR)*(1.03+1.35*(Xt(V_index))-Xt(R_index)+(Xt(V_index)+0.38)^2*3.3*PIRoff);
end

function xfprime = Fderivative(t,Xt,F_index, op_V_index)
  global synThres TauSyn;
  % Hstep is determined by whether the opposing ventral level 
  Hstep = (Xt(op_V_index) > synThres);
  xfprime=(1/TauSyn)*(-Xt(F_index) + Hstep);
end

function xgprime = Gderivative(t,Xt,G_index, F_index)
  global TauSyn;
  xgprime = (1/TauSyn)*(-Xt(G_index) + Xt(F_index));
end

function V=Spikes(Vin)
  Last = size(Vin)(2);
  V = [0 (Vin(1:Last - 1) < -0.2).*(Vin(2:Last) >= -0.2)];
end

%return 0 for system not giving out of phase in sync spiking
function Hz=analyzeSeries(XT)
  Time = 1:size(XT)(2);
  spikesD = Spikes(XT(2,:))
  spikesV = Spikes(XT(3,:));
  nSpikesD = sum(spikesD);
  nSpikesV = sum(spikesV);
  if(nSpikesD<3 || nSpikesV < 3 || abs(nSpikesD - nSpikesV) > 2)
    %CPG not working
    Hz = 0;
  end
  
  plot(Time,spikesD,Time,spikesV)
  %count spikes blue
  %count spikes red
  %verify that they are out of phase and the number of spikes not more than 2 off

end

