% fwd Euler for Clione swimming

function XT=Clione (kin,dorsalStimin,noiseLevel,doPlot)
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
  if(nargin<4)
    doPlot = 1;
  end
  if(nargin<3)
    noiseLevel = 0;
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

  %Noise amplitude is average 0.068788 without scaling
  avgNoise = 0.068788;
  if noiseLevel != 0
    Noise = [];
    for i = 1:8
      Noise = [Noise; pinkNoise(Last)];
    end
    Noise = Noise.*noiseLevel;
  else
    Noise = zeros(8,Last);
  end

  % Forward Euler:
  for t = 2:Last
    X(:,t) = X(:,t-1) + (derivative(X(:,t-1),t-1) + Noise(:,t-1)).*DT;
  end;

  %The next two lines use Octave's built-in solver.
  %[X,istate,msg] = lsode(@derivative,Xinit,Time);
  %X = X';

  if(doPlot)
    figure(1), ZA = plot(Time, 100*X(1,:), 'r', Time, 100*X(2,:).-150, 'b-'); set(ZA, 'LineWidth', 2);
    ylabel('V (mV'); xlabel('Time (ms)');
    title(cstrcat("Clione at k=",num2str(k),", stimulus of I=",num2str(dorsalStim),
                  ", noise level of ", num2str(noiseLevel*avgNoise)," "));

    VV = -0.9:0.01:1.5;
    DVdt = -0.5*((1.37 + 3.67*VV + 2.51*VV.^2).*(VV - 0.55) - dorsalStim/13)./(VV + 0.92);
    DRdt = 1.35*VV + 1.03;
    figure(2), ZB = plot(VV, DVdt, 'k-', VV, DRdt, 'b-', X(1,:), X(3,:), 'r-'); axis([-1, 0.6, 0, 1]);
    set(ZB, 'LineWidth', 2); axis square;
  end
  XT = [Time; X];
end


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

function MHz=SpikeRate(SpkTime)
  Final = length(SpkTime);
  Rates = SpkTime(2:Final) - SpkTime(1:Final - 1);
  Leng = length(Rates);
  MHz= mean(Rates(floor(Leng/2):Leng));
end

%return 0 for system not giving out of phase in sync spiking
function MHz=analyzeSeries(XT)
  Time = 1:size(XT)(2);
  SpikesD = Spikes(XT(2,:));
  SpikesV = Spikes(XT(3,:));
  nSpikesD = sum(SpikesD);
  nSpikesV = sum(SpikesV);
  if(nSpikesD<3 || nSpikesV < 3 || abs(nSpikesD - nSpikesV) > 2)
    %CPG not working
    MHz = 0;
    return
  end
  SpikeTimesD = zeros(1, nSpikesD);
  SpikeTimesV = zeros(1, nSpikesV);
  nD = 1; nV = 1; DT = XT(1,2)-XT(1,1);
  for T = 1:length(SpikesD);  %record times of spikes in spike time arrays
  	if SpikesD(T) == 1; 
      SpikeTimesD(nD) = T*DT;
      nD = nD + 1; 
    end;
  	if SpikesV(T) == 1; 
      SpikeTimesV(nV) = T*DT; nV = nV + 1; 
    end;
  end;
  SpikeRateD = SpikeRate(SpikeTimesD);
  SpikeRateV = SpikeRate(SpikeTimesV);
  %verify that they are out of phase
  flag = true;
  try
    for i = 3:nSpikesD
      flag = flag && SpikeTimesD(i-1) + SpikeRateD/2 - SpikeTimesV(i-1) <= 10*DT;
      flag = flag && SpikeTimesV(i-1) + SpikeRateD/2 - SpikeTimesD(i)<=10*DT;
    end
  catch
    % nSpikesD and nSpikesV are too different
    flag = false;
  end

  if flag
    MHz = (SpikeRateD + SpikeRateV)/2;
  else
    MHz = 0;
  end

  %plot(Time,SpikesD,Time,SpikesV)
end

function [kRange,noiseRange,HeatMap] = gridSearch()
  kRange = [0:10 12:2:30 35:5:600];
  Stim = 0.5;
  noiseRange =  0:0.002:0.01;
  HeatMap = zeros(length(kRange),length(noiseRange));
  for kInd = 1:length(kRange)
    k = kRange(kInd);
    for noiseInd = 1:length(noiseRange)
      noise = noiseRange(noiseInd);
      HeatMap(kInd,noiseInd) = analyzeSeries(Clione(k,Stim,noise,false));
      disp(noiseRange(noiseInd));
    end
    disp(k);
  end
  save noiseData kRange noiseRange HeatMap
end

function x=pinkNoise(Nx)
  % Nx number of samples to synthesize
  B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
  A = [1 -2.494956002   2.017265875  -0.522189400];
  nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
  v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
  x = filter(B,A,v);    % Apply 1/F roll-off to PSD
  x = x(nT60+1:end);    % Skip transient response
end

function graphHeatMap(noiseRange,kRange,HeatMap)
  avgNoise = 0.068788;
  contourf(noiseRange.*avgNoise, kRange, HeatMap,[15,14,13,12,11,10,9,0]);
  colorbar
  ylabel('k (connection strength)');
  xlabel('Noise level: mean mV');
end
%
%function HeatColumn = stimSearch(k)
%  stimRange = 0:0.1:2;
%  HeatColumn = zeros(1,length(stimRange));
%  for stimInd = 1:length(stimRange)
%    HeatColumn(stimInd) = analyzeSeries(Clione(k,stimRange(stimInd),false));
%  end
%end
%
%function [kRange,stimRange,HeatMap] = pGridSearch(nProcs)
%  % k is always the row number + 1
%  kRange = 0:600;
%  stimRange = 0:0.5:2;
%  HeatMap = zeros(length(kRange),length(stimRange));
%  for k = kRange
%    HeatMap(k+1,:) = stimSearch(k);
%  end
%end
