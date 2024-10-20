function networkAnalysisSkipSpk
clc;clear all;close all
global NoConn neuronModel

format long

NoConn=-1;
neuronModel=1;

tBegin=0;tEnd=100;dt=0.0001;

tol_eq=eps;
rand('seed',5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The number of neurons in the network.
N=6; % Osc. is ordered from 1 to 100.

% Excitatory index.
N_ex=4; % A neuron is excitatory iff its index is greater or equal N_ex.

% The period of pattern that is normalized by the gamma period (30 ms).
Gamma_period=30; % [ms]

% Parameters of LIF oscillator.
gamma=1;
% gamma=-1.5;
% I=1.2;
I=2;

dataMode=2;

% Period of each oscillator.
T=1;

% Type of oscillator: -1 (inhibitory) and 1 (excitatory).
osc_types=ones(1,N);
osc_types(1,1:1:(N_ex-1))=-1;

% Upper and lower bound of weights.
lbInh=-Inf;
ubInh=-0.0001;
lbExc=0.01;
ubExc=Inf;

% tol_boundary
tol_boundary=0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Natural_period=(-1/gamma)*log((1-I/gamma)/(-I/gamma));
display(strcat('Natural period of osc.:',num2str(Natural_period),'=',num2str(Natural_period.*Gamma_period),' ms'));

periods=ones(1,N).*T;
gamma_i=ones(1,N).*gamma;
I_i=ones(1,N).*I;

% Generate connection matrix from osc. j (column) to osc. i (row).
conns=genConnMatrix(N,0);

% Spike time of each oscillator.
spk_times=genSpikeTime(N,N_ex,1,1); % [degree] from 0 to 360.
spk_times=spk_times.*Gamma_period./360; % [ms]
spk_times=spk_times/Gamma_period; %[-] normalized by the gamma period

center_Spk_Ex=0.158237860621827;
center_Spk_Inh=0.245569933069455+0.2;
% spk_times=[0.33 0.3461 0.35 0.13 0.1354 0.14];
% spk_times=[0.24 0.3461 0.44 0.03 0.1354 0.23];

i=0;
for dInh=0.01:0.001:0.10
    i=i+1;
    j=0;
    for dExc=0.01:0.001:0.10
        j=j+1;    
    
spk_times=[center_Spk_Inh-dInh center_Spk_Inh center_Spk_Inh+dInh center_Spk_Ex-dExc center_Spk_Ex center_Spk_Ex+dExc];

% Show spike times 
figure(1)
statSpk_Times(spk_times,N_ex,N,Gamma_period);
ylim([0 N+1]);

spk_times_ref=...
    [0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827
    0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827
    0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827
    0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827
    0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827
    0.435569933069455   0.445569933069455   0.545569933069455   0.058237860621827   0.158237860621827   0.258237860621827];


conns = ...
     [1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     0     0     0
     1     1     1     1     0     0
     1     1     1     1     1     0];


delays_ref = ...
   [0.360126138366818   0.595890195139633   0.844092978795102   0.377332072447628   0.277332072447628   0.177332072447628
   0.650007531102750   0.289403749298958   0.071488930737362   0.387332072447628   0.287332072447628   0.187332072447628
   0.716575243923150   0.108814467586956   0.016065447598726   0.487332072447628   0.387332072447628   0.287332072447628
   0.442171823182881   0.088687808202446   0.006140629694541                 NaN                 NaN                 NaN
   0.046559640565305   0.660235439448370   0.190384502558740   0.100000000000000                 NaN                 NaN
   0.684370528147573   0.384741880250650   0.512301905662998   0.200000000000000   0.100000000000000                 NaN];

abs_Arrive_Time=spk_times_ref+delays_ref;
rel_Arrive_Time=abs_Arrive_Time-spk_times_ref';

spk_times_Mat(1,:)=spk_times;
spk_times_Mat(2,:)=spk_times;
spk_times_Mat(3,:)=spk_times;
spk_times_Mat(4,:)=spk_times;
spk_times_Mat(5,:)=spk_times;
spk_times_Mat(6,:)=spk_times;

delays=spk_times_Mat'+rel_Arrive_Time-spk_times_Mat;

lbInh=conns.*lbInh;
ubInh=conns.*ubInh;
lbExc=conns.*lbExc;
ubExc=conns.*ubExc;

% Get weights.
display('Calculating weights');
[weights,arr4eachOsc,v0,ret_msg]=getWeightsV3(size(spk_times,2),spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary);
display(ret_msg);

figure(2)
calCheckCalculatedWeightsV3(weights,arr4eachOsc,periods,gamma_i,I_i,spk_times,T);

if (strcmp(ret_msg,'Succes!!!')==1)    
%   Simulate network.
    figure(3)
    display('Simulating network');    
    tol_eq=1e-10;
%     simNetworkEBSkipSpk(weights, delays, v0, spk_times, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, tol_eq, osc_types);    
    d_Ref_Pert(i,j)=simNetworkEBSkipSpkShort(weights, delays, v0, spk_times, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, tol_eq, osc_types)
end

maximize('all');

    end
end

save('networkAnalysisSkipSpk.mat',d_Ref_Pert);

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

function statSpk_Times(spk_times,N_ex,N,Gamma_period)

hold on
plot(spk_times(1,N_ex:1:N).*Gamma_period,N_ex:1:N,'*r');
plot(spk_times(1,1:1:(N_ex-1)).*Gamma_period,1:1:(N_ex-1),'*b');
xlim([0 Gamma_period]);
xlabel('time [ms]');
ylabel('oscillator');
grid on

end

function statDelays(delays,N_ex,Gamma_period)
E_E=0;E_E_n=0;
E_I=0;E_I_n=0;
I_E=0;I_E_n=0;
I_I=0;I_I_n=0;

N=size(delays,2);

for i=1:1:N
    for j=1:1:N
        if ~isnan(delays(i,j))
            if ((i>=N_ex) && (j>=N_ex))
                % E to E
                E_E_n=E_E_n+1;
                E_E(E_E_n)=delays(i,j);                
            elseif ((i>=N_ex) && (j<N_ex))
                % I to E
                I_E_n=I_E_n+1;
                I_E(I_E_n)=delays(i,j);                
            elseif ((i<N_ex) && (j>=N_ex))
                % E to I
                E_I_n=E_I_n+1;
                E_I(E_I_n)=delays(i,j);                
            elseif ((i<N_ex) && (j<N_ex))
                % I to I
                I_I_n=I_I_n+1;
                I_I(I_I_n)=delays(i,j);                
            end
        end
    end
end

display(strcat('E->E: mean:',num2str(mean(E_E)*Gamma_period),' ms std:',num2str(std(E_E)*Gamma_period),' ms'));
display(strcat('I->E: mean:',num2str(mean(I_E)*Gamma_period),' ms std:',num2str(std(I_E)*Gamma_period),' ms'));
display(strcat('E->I: mean:',num2str(mean(E_I)*Gamma_period),' ms std:',num2str(std(E_I)*Gamma_period),' ms'));
display(strcat('I->I: mean:',num2str(mean(I_I)*Gamma_period),' ms std:',num2str(std(I_I)*Gamma_period),' ms'));
end