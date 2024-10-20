function calGetWeights
clc;clear all;close all
format long

%% Spike times of the oscillators.
spk_times=[0, 0.0751724, 0.106362, 0.484572, 0.726624, 1.28886];
% spk_times=conv2TrueTimeLIF(spk_times);

%% Parameters of each LIF oscillator.
gamma_i=[1,1,1,1,1,1];
I_i=[1.2,1.2,1.2,1.2,1.2,1.2];

%% Parameters of the connections.
% 1) Homogenous delays
conns=[ 1	0	1	1	1	1;
        1	1	1	1	1	1;
        0	1	1	1   1   1;
        1   1   0	1   0	0;
        1   0	1   1   1   0;
        0	1   1   1   0	1];
delay=0.125;    
conns(conns==0)=NaN;
delays=conns*delay;
% delays= ...
% [
%    0.125000000000000                 NaN   0.125000000000000   0.725000000000000   0.625000000000000   0.125000000000000
%    0.125000000000000   0.125000000000000   0.125000000000000   0.825000000000000   0.625000000000000   0.125000000000000
%                  NaN   0.125000000000000   0.125000000000000   0.725000000000000   0.625000000000000   0.125000000000000
%    0.625000000000000   0.625000000000000                 NaN   1.225000000000000                 NaN                 NaN
%    0.925000000000000                 NaN   0.725000000000000   0.125000000000000   1.225000000000000                 NaN
%                  NaN   0.125000000000000   0.125000000000000   0.625000000000000                 NaN   1.225000000000000
% ];
% delays=conv2TrueTimeLIF(delays);

%% Period of the pattern.
T=1.5;
% T=conv2TrueTimeLIF(T);

%% Osc. is excitatory iff its index is greater than or equal ex_neuron_i.
ex_neuron_i=4;

%% Upper and lower bound of weights.
lbInh=-Inf;
ubInh=-0.01;
lbExc=0.01;
ubExc=Inf;

%% tol_boundary
tol_boundary=0.01;

%% Get weights.
display(strcat('Spike times:',num2str(spk_times)));
display(strcat('Period:',num2str(T)));
[weights,arr4eachOsc,v0,ret_msg]=getWeights(spk_times, gamma_i, I_i, delays, T, ex_neuron_i,0.000001,lbInh,ubInh,lbExc,ubExc,tol_boundary)

calCheckCalculatedWeightsV2(weights,arr4eachOsc,T,gamma_i(1,1),I_i(1,1));

end

function true_time=conv2TrueTimeLIF(dimensionless_time)
    I=1.2;
    T=log(I/(I-1));
    
    true_time=dimensionless_time.*T;
end

