function calSimNetwork
clc
clear all
close all

global NoConn
NoConn=-1;

tEnd=10;
dt=0.0001;

%% Data from Raoul
% tSpikes=[0, 0.0751724, 0.106362, 0.484572, 0.726624, 1.28886];
% tau=0.125;  % Dimensionless time delay between the oscillators. It's scaled by the natural period of the oscillator.
% T=1.5;      % Dimensionless Period of the pattern. It's scaled by the natural period of the oscillator.
% % It's the weight from osc. identified by column to osc. identified by row.
% epsiMat=[
%     0.224596	0	-0.300411	-0.686132	-0.241245	0.126174;
%     -0.902917	1.10085	-1.17011	-1.20955	-0.757655	0.51012;
%     0	0.672664	-1.1734	-0.879224	-0.266808	0.202992;
%     -0.0806747	-0.0726107	0	-0.380294	0	0;
%     -1.25423	0	-1.59692	1.20045	-0.0867559	0;
%     0	-0.878748	0.445186	-0.380777	0	0.615526
%     ];
% % Phases of the oscillators just after arrival of the spike from osc. 1.
% initVoltages=phi2vLIF([0.273885, -0.285216, 0.0186377, 0.777742, -0.130416, 0.908975]);
% tauMat=conv2TrueTimeLIF((epsiMat./epsiMat).*tau);    
% dataMode=0;
% tBegin=conv2TrueTimeLIF(tau);
% tSpikes=conv2TrueTimeLIF(tSpikes);
% T=conv2TrueTimeLIF(T);

%% Data1 from Atthaphon
% tSpikes=[0, 0.0751724, 0.106362, 0.484572, 0.726624, 1.28886];
% tau=0.125;  % Dimensionless time delay between the oscillators. It's scaled by the natural period of the oscillator.
% T=1.5;      % Dimensionless Period of the pattern. It's scaled by the natural period of the oscillator.
% epsiMat=[
%   -0.010000000000000                 NaN  -0.010000000000000   0.025495450742691   0.029386989253662   0.031274029791004
%   -0.010000000000000  -0.010000000000000  -0.010000000000000   0.029559293509389   0.030828827054042   0.032808452998963
%                  NaN  -0.010000000000000  -0.010000000000000   0.028216467762187   0.032523333016240   0.034611769068037
%   -0.010000000000000  -0.010000000000000                 NaN   0.096227947025910                 NaN                 NaN
%   -0.010000000000000                 NaN  -0.010000000000000   0.047430723882196   0.040500819233264                 NaN
%                  NaN  -0.010000000000000  -0.010000000000000   0.048879907383227                 NaN   0.044418427631647
%     ];
% epsiMat(isnan(epsiMat))=0;
% initVoltages=[0   0.984402042687980   0.977569902502031   0.756649760867523   0.635484195382876   0.228378543974360];
% % tauMat=(epsiMat./epsiMat).*tau; 
% tauMat= ...
% [
%    0.125000000000000                 NaN   0.125000000000000   0.725000000000000   0.625000000000000   0.125000000000000
%    0.125000000000000   0.125000000000000   0.125000000000000   0.825000000000000   0.625000000000000   0.125000000000000
%                  NaN   0.125000000000000   0.125000000000000   0.725000000000000   0.625000000000000   0.125000000000000
%    0.625000000000000   0.625000000000000                 NaN   1.225000000000000                 NaN                 NaN
%    0.925000000000000                 NaN   0.725000000000000   0.125000000000000   1.225000000000000                 NaN
%                  NaN   0.125000000000000   0.125000000000000   0.625000000000000                 NaN   1.225000000000000
% ];
% dataMode=1;
% tBegin=0;
% tSpikes=conv2TrueTimeLIF(tSpikes);
% T=conv2TrueTimeLIF(T);

%% Data2 from Atthaphon
% tSpikes=[0, 0.0751724, 0.106362, 0.484572, 0.726624, 1.28886];
% tau=0.125;  % Dimensionless time delay between the oscillators. It's scaled by the natural period of the oscillator.
% T=1.5;      % Dimensionless Period of the pattern. It's scaled by the natural period of the oscillator.
% epsiMat=[
%   -0.049051384589803                 NaN  -0.059349490384739  -0.116874987804962  -0.180332766390298  -0.032584048341170
%   -0.040283574292665  -0.046091781545172  -0.048740919856572  -0.095983880854084  -0.148098742839474  -0.066118717995315
%                  NaN  -0.044931763413204  -0.047514229351085  -0.093568199827343  -0.144371457383019  -0.081413967334527
%   -0.077786905679758  -0.089002456376348                 NaN  -0.282421053264095                 NaN                 NaN
%   -0.094436270134724                 NaN  -0.114262717621566  -0.045829578776834  -0.023623014927433                 NaN
%                  NaN  -0.122173246687944  -0.129195188942702  -0.254419811092779                 NaN  -0.073144514753341
%     ];
% epsiMat(isnan(epsiMat))=0;
% initVoltages=[0   0.971159097023676   0.958003403081906   0.948191950605179   0.892449956526363   0.315300141034882];
% tauMat=conv2TrueTimeLIF((epsiMat./epsiMat).*tau);    
% dataMode=1;
% tBegin=0;
% tSpikes=conv2TrueTimeLIF(tSpikes);
% T=conv2TrueTimeLIF(T);

%% Data3 from Atthaphon
tSpikes=[ 
    0.243885563697163   0.244110779026641   0.244220325282536   0.169023264095047   0.187099903977570   0.181702733959407
    ];
T=1;      % Dimensionless Period of the pattern. It's scaled by the natural period of the oscillator.
epsiMat=[
  -0.270813197023367  -0.298336916067365  -0.280587960493016   0.100000000000000   0.100000000000000   0.100000000000000
  -0.368661739349888  -0.214785104723580  -0.241220540515445   0.100000000000000   0.100000000000000   0.100000000000000
  -0.263653523689332  -0.350644659066184  -0.254278954991713   0.100000000000000   0.100000000000000   0.100000000000000
  -0.123286058115782  -0.144822996901053  -0.199616297799992                 NaN                 NaN                 NaN
  -0.212750973897415  -0.201976993372070  -0.233491894562677   0.100000000000000                 NaN   0.100000000000000
  -0.192502954479574  -0.244903072309463  -0.184554602266305   0.100000000000000                 NaN                 NaN
  ];
initVoltages=[0.340912606838216   0.747962161080066   0.895392921638137   0.815777302956668   0.553109463691007   0.680723106365669];
tauMat=[
   0.544050455346727   0.640619405396694   0.579173680892507   0.074862299602115   0.056785659719593   0.062182829737756
   0.856228227770579   0.315761499719583   0.431725559005492   0.075087514931594   0.057010875049072   0.062408045067235
   0.842964859154633   0.003635333290677   0.806426179039491   0.075197061187489   0.057120421304966   0.062517591323129
   0.077757681756500   0.238537421233397   0.559312602107313                 NaN                 NaN                 NaN
   0.590236156218180   0.538042392501613   0.682926413241400   0.018076639882522                 NaN   0.005397170018163
   0.312064885780965   0.552590677934119   0.269563989579371   0.012679469864360                 NaN                 NaN 
    ];    
dataMode=2;
tBegin=0;
gamma=1;
I=2;
gamma_i=ones(1,size(tSpikes,2)).*gamma;
I_i=ones(1,size(tSpikes,2)).*I;
osc_types=ones(1,size(tSpikes,2));

%% Begin processing
dV=rand(1,size(tSpikes,2))*0; % perturb the initial voltages.
% simNetwork(epsiMat, tauMat, initVoltages, tSpikes, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, eps);
simNetworkEB(epsiMat, tauMat, initVoltages, tSpikes, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, eps, osc_types, 0);

%% Draw reference lines
noOscs=size(tSpikes,2);
plotGrid(T,noOscs,tEnd);

maximize('all');
end

function plotGrid(T,noOscs,tEnd)
    y=linspace(0,noOscs+1,noOscs*100);
    T_truetime=T;
    i=1;
    while (i*T_truetime<=tEnd)
        figure(1)
        hold on
        plot(i*T_truetime,y,'r')
        i=i+1;
    end        
end

function true_time=conv2TrueTimeLIF(dimensionless_time)
    I=1.2;
    T=log(I/(I-1));
    
    true_time=dimensionless_time.*T;
end

function voltage=phi2vLIF(phi)
    gamma=1;
    I=1.2;
    T=log(I/(I-1));
    
    voltage=(I./gamma).*(1-exp(-gamma.*phi.*T));
end

function phi=v2phiLIF(v)
    gamma=1;
    I=1.2;
    T=log(I/(I-1));
    
    phi=-(1/(gamma*T)).*log(1-gamma*v/I);
end