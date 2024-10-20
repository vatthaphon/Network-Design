function calGetWeightsV2
clc;clear all;close all
format long

%% Parameters of each LIF oscillator.
gamma_i=[1,1,1,1,1,1];
I_i=[1.2,1.2,1.2,1.2,1.2,1.2];
Natural_period=log(I_i./(I_i-1));

%% Spike times of the oscillators.
spk_times=[0, 0.0751724, 0.106362, 0.484572, 0.726624, 1.28886];
spk_times=rand(1,6);

%% Parameters of the connections.
% 1) Homogenous delays
% conns=[ 1	0	1	1	1	1;
%         1	1	1	1	1	1;
%         0	1	1	1   1   1;
%         1   1   0	1   0	0;
%         1   0	1   1   1   0;
%         0	1   1   1   0	1];
conns=[ 1	1	1	1	1	1;
        1	1	1	1	1	1;
        1	1	1	1   1   1;
        1   1   1	1   1	1;
        1   1	1   1   1   1;
        1	1   1   1   1	1];
conns(conns==0)=NaN;
delay=0.125;    
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

% delays=rand(6,6);

%% Period of each oscillator.
% T=1.5;
T=1;
% periods=[T T T T T T];
periods=[T T T T T T];

%% Type of oscillator: -1 (inhibitory) and 1 (excitatory).
osc_types=[-1 -1 -1 1 1 1];

%% Upper and lower bound of weights.
lbInh=conns.*(-Inf);
ubInh=conns.*(-0.001);
lbExc=conns.*0.001;
ubExc=conns.*Inf;

%% tol_boundary
tol_boundary=0.01;

%% Get weights.
display(strcat('Spike times:',num2str(spk_times)));
display(strcat('Period:',num2str(T)));
[weights,arr4eachOsc,v0,ret_msg]=getWeightsV2(size(spk_times,2),spk_times,periods,gamma_i,I_i,osc_types,delays,eps,lbInh,ubInh,lbExc,ubExc,tol_boundary)

calCheckCalculatedWeightsV3(weights,arr4eachOsc,periods,gamma_i,I_i,spk_times,T);

end

