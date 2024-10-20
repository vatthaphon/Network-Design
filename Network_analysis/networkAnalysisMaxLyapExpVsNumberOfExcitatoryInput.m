function networkAnalysisMaxLyapExpVsNumberOfExcitatoryInput
clc;clear all;close all

addpath('../Conn_matrix');
addpath('../Check_calculated_weights');
addpath('../Extract_pdf_from_the_paper');
addpath('../Gen_delays');
addpath('../Get_weights');
addpath('../LyapExp');
addpath('../Sim_network');
addpath('../Network_analysis');

global NoConn neuronModel

format long

eclipse_ws='C:\Attha\Radboud U\Source code\Attha\paper2_Raoul\eclipse_ws\';
% eclipse_ws='C:\Users\attha\Documents\Work\Radboud U\Source code\Attha\paper2_Raoul\eclipse_ws\';

NoConn=-1;
neuronModel=1;

tBegin=0;tEnd=2;dt=0.0001;

tol_eq=eps;

for i=10:1:49
% for i=50:-1:11
    Rnd_Number=6;
    rand('seed',Rnd_Number);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileName=strcat('Gamma_notsatisfycond_ExcInputLeft_',num2str(i),'_50');
%     fileName=strcat('Gamma_notsatisfycond_ExcInputLeft_10_kill_',num2str(i));
    
    % The number of neurons in the network.
    N=50; % Osc. is ordered from 1 to 100.
    
    % Excitatory index.
    N_ex=10; % A neuron is excitatory iff its index is greater or equal N_ex.
    
    % The period of pattern that is normalized by the gamma period (30 ms).
    Gamma_period=30; % [ms]
    
    % Parameters of LIF oscillator.
    gamma=1;
%     I=2.0;
    I=1.2;
    
    dataMode=2;
    
    % Period of each oscillator.
    T=1;
    
    % Type of oscillator: -1 (inhibitory) and 1 (excitatory).
    osc_types=ones(1,N);osc_types(1,1:1:(N_ex-1))=-1;
    
    % Upper and lower bound of weights.
    lbInh=-Inf;
    ubInh=-0.0001;
    lbExc=0.01;
    ubExc=Inf;
    
    ubInh=-0.01;
    lbExc=0.01;
    
    % ubInh=-0.0001;
    % lbExc=0.0001;
    
    % Prohibit accidentally spiking before we want.
    tol_boundary=0.05;
    
    % To make sure that the excitatory inputs bring the state more than one.
    % This option is only used when the suddenly spking mechanism.
    tol_spk=0.0;
    
    % The inh. input should arrive at phase<tMax: the exc. input arrive at
    % phase>=tMax.
    tMax=0.6;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Natural_period=(-1/gamma)*log((1-I/gamma)/(-I/gamma));
    display(strcat('Natural period of osc.:',num2str(Natural_period),'=',num2str(Natural_period.*Gamma_period),' ms'));
    
    periods=ones(1,N).*T;
    gamma_i=ones(1,N).*gamma;
    I_i=ones(1,N).*I;
    
    % Spike time of each oscillator.
    spk_times=genSpikeTime(N,N_ex,1,1); % [degree] from 0 to 360.
%     spk_times=genSpikeTime(N,N_ex,10000,10000000); % [degree] from 0 to 360.
    spk_times=spk_times.*Gamma_period./360; % [ms]
    spk_times=spk_times/Gamma_period; %[-] normalized by the gamma period
    
%     % Show spike times
%     figure(1)
%     statSpk_Times(spk_times,N_ex,N,Gamma_period);
%     ylim([0 N+1]);
    
    %% Generate connection matrix from osc. j (column) to osc. i (row).
    conns=genConnMatrix(N,0);
        
    %% Generate delay for postsynaptic neurons suddenly spike
    % Not allow unreasonable connections.
    % % Get rid of unreasonable connections, i.e., spike time of presynaptic
    % % excitatory neuron should be earier than spike time of postsynaptic neuron.
    % conns=killConns(conns,spk_times,N_ex);
    %
    % % Determine delay times for each connection
    % [conns,delays]=genDelays(conns,spk_times,N_ex,1,tol_eq);
    
    % % Allow unreasonable connections.
    % delays=genDelaysV1(conns,spk_times,N_ex,1,tol_eq);
    
    %% Generate delay for Postsynaptic neurons donot spike suddenly
    % Allow unreasonable connections.
    delays=genDelaysV2(conns,spk_times,N_ex,1,tol_eq,tMax);
    
    % % Get delays that provide the same dynamics as reference
%     delays=getDelayWithSameRelativeArr(N,spk_times,T,tol_eq,arr4eachOsc);

%     conns(:,10:1:i)=0;
%     delays(:,10:1:i)=NaN;

    conns(:,i:1:50)=0;
    delays(:,i:1:50)=NaN;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Statistic information for delay
    statDelays(delays,N_ex,Gamma_period);
    
    lbInh=conns.*lbInh;
    ubInh=conns.*ubInh;
    lbExc=conns.*lbExc;
    ubExc=conns.*ubExc;
    
    % Get weights.
    display('Calculating weights');
    
    % Postsynaptic neurons suddenly spike.
    [weights,arr4eachOsc,v0,ret_msg]=getWeightsV4(size(spk_times,2),spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary,tol_spk);
    
    if (isWeightCorrect(weights,osc_types,N,tol_eq,ret_msg)==0)
        display('Matlab doesnot provide correct weights');
    end
    
    display(ret_msg);
    
    % figure(2)
    % calCheckCalculatedWeightsV3(weights,arr4eachOsc,periods,gamma_i,I_i,spk_times,T);
        
    % Save calculated data
%     save(strcat('../Data/Gamma_notsatisfycond_ExcInputLeft/',fileName,'.mat'),'N','tol_spk','lbInh','ubInh','lbExc','ubExc','conns','tol_boundary', ...
%         'Gamma_period','N_ex','fileName','Rnd_Number','NoConn', ...
%         'neuronModel','weights','arr4eachOsc','delays','v0', ...
%         'spk_times','tBegin','tEnd','dt','dataMode','T','periods', ...
%         'gamma_i','I_i','tol_eq','osc_types');
%     mysave(eclipse_ws,fileName,N,NoConn,weights,delays,v0,spk_times,T,tol_eq,gamma_i,I_i);
    
% if (strcmp(ret_msg,'Succes!!!')==1)
%     figure(3)
%     display('Simulating network');
%     dV=rand(1,N)*0; % perturb the initial voltages.
%     tol_eq=1e-14;
%     simNetworkEB(weights, delays, v0, spk_times, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, tol_eq, osc_types, dV);    
% end    
% 
% k = waitforbuttonpress;
% 
% clf;

end

end

function mysave(eclipse_ws,fileName,N,NoConn,weights,delays,v0,spk_times,T,tol_eq,gamma_i,I_i)
    %Replace data
    delays(isnan(delays)) = NoConn;
    delays(isnan(weights)) = NoConn;
            
    fid = fopen(strcat(eclipse_ws,'Lyap_Exp\data\',fileName),'w');                    
    
    fwrite(fid, N, 'int');
    fwrite(fid, NoConn, 'int');    
    fwrite(fid, weights, 'double');
    fwrite(fid, delays, 'double');
    fwrite(fid, v0, 'double');
    fwrite(fid, spk_times, 'double');
    fwrite(fid, T, 'double');
    fwrite(fid, tol_eq, 'double');
    fwrite(fid, gamma_i, 'double');
    fwrite(fid, I_i, 'double');
    
    fclose(fid);
end

function val=isWeightCorrect(weights,osc_types,N,tol_eq,ret_msg)
val=1;
for i=1:1:N
    for j=1:1:N
        tmp=weights(i,j)*osc_types(1,j);
        if ~isnan(tmp)
            if (cmp(tmp,0,tol_eq)<0)
                val=0;
                return;
            end
        end
    end
end

if (strcmp(ret_msg,'Succes!!!')==0)
    val=0;
    return
end

end

function delays=getDelayWithSameRelativeArr(nOscs,spk_times,T,tol_eq,arr4eachOsc)
delays=NaN(nOscs);
for i=1:1:nOscs
    for j=1:1:nOscs
        if ~isnan(arr4eachOsc(i,j))
            if (cmp(spk_times(1,j),spk_times(1,i),tol_eq)<0)
                delays(i,j)=arr4eachOsc(i,j)+spk_times(1,i)-spk_times(1,j);
            else
                if (cmp(arr4eachOsc(i,j)+spk_times(1,i),spk_times(1,j),tol_eq)<=0)
                    delays(i,j)=arr4eachOsc(i,j)+spk_times(1,i)+T-spk_times(1,j);
                else
                    delays(i,j)=arr4eachOsc(i,j)+spk_times(1,i)-spk_times(1,j);
                end
            end
        end
    end
end

end

function val=cmp(x, y, tol_eq)
% CMP Two-value comparison
%   val = cmp(x, y, tol_eq)
% Input
%   x           the first number.
%   y           the second number.
%   tol_eq      if the first and second numbers are different less than
%               tol_eq, we say that the two numbers are equal.
% Output
%   val         0   : two numbers are the same.
%               -1  : the first number is less than the second number.
%               1   : the first number is greater than the second number.

if  (abs(x-y)<tol_eq)
    val=0;
elseif (x<y)
    val=-1;
else
    val=1;
end

end

function [n,r]=mMod(spk_time,T,x,tol_eq)
if isnan(x)
    n=NaN;r=NaN;
else
    if (cmp(x,spk_time,tol_eq)<0)
        n=0;
        r=T-(spk_time-x);
    else
        n=0;
        while (1)
            if (cmp(x,spk_time+(n+1)*T,tol_eq)<0)
                r=x-(spk_time+n*T);
                
                % If an input make a post neuron spike suddenly, we say
                % that the input arrives at T instead of 0.
                if (cmp(r,0,tol_eq)==0)
                    r=T;
                end
                return
            end
            n=n+1;
        end
    end
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