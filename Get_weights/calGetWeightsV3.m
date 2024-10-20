function calGetWeightsV3
clc;clear all;close all
format long

tol_eq=eps;
rand('seed',2);

%% The number of neurons in the network.
N=6; % Osc. is ordered from 1 to 100.

%% Excitatory index.
N_ex=4; % A neuron is excitatory iff its index is greater or equal N_ex.

%% Generate connection matrix from osc. j (column) to osc. i (row).
conns=genConnMatrix(N,0);

%% The period of pattern that is normalized by the gamma period (30 ms).
Gamma_period=30; % [ms]

%% Spike time of each oscillator.
spk_times=genSpikeTime(N,N_ex,1,1); % [degree] from 0 to 360.
spk_times=spk_times.*Gamma_period./360; % [ms]
spk_times=spk_times/Gamma_period %[-] normalized by the gamma period

% Show spike times 
statSpk_Times(spk_times,N_ex,N,Gamma_period);

%% Get rid of unreasonable connections, i.e., spike time of presynaptic 
% excitatory neuron should be earier than spike time of postsynaptic neuron.
conns=killConns(conns,spk_times,N_ex);

%% Determine delay times for each connection
[conns,delays]=genDelays(conns,spk_times,N_ex,1,tol_eq)

% Statistic information for delay
statDelays(delays,N_ex,Gamma_period);

%% Parameters of each LIF oscillator.
gamma=1;
% I=1.2;
I=2;
gamma_i=ones(1,N).*gamma;
I_i=ones(1,N).*I;
Natural_period=log(I./(I-1));
display(strcat('Natural period of osc.:',num2str(Natural_period),'=',num2str(Natural_period.*Gamma_period),' ms'));

%% Period of each oscillator.
T=1;
periods=ones(1,N).*T;

%% Type of oscillator: -1 (inhibitory) and 1 (excitatory).
osc_types=ones(1,N);
osc_types(1,1:1:(N_ex-1))=-1;

%% Upper and lower bound of weights.
lbInh=conns.*(-Inf);
ubInh=conns.*(-0.1);
lbExc=conns.*0.1;
ubExc=conns.*Inf;

%% tol_boundary
tol_boundary=0.05;

%% Get weights.
[weights,arr4eachOsc,v0,ret_msg]=getWeightsV3(size(spk_times,2),spk_times,periods,gamma_i,I_i,osc_types,delays,eps,lbInh,ubInh,lbExc,ubExc,tol_boundary)

calCheckCalculatedWeightsV3(weights,arr4eachOsc,periods,gamma_i,I_i,spk_times,T);


end

function statSpk_Times(spk_times,N_ex,N,Gamma_period)

figure(1);hold on
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