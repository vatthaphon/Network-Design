function networkAnalysis
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

tBegin=0;tEnd=1;dt=0.0001;

tol_eq=eps;
Rnd_Number=6;
% Rnd_Number=10;
% Rnd_Number=20;
rand('seed',Rnd_Number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filename
% fileName='NotPureInh';
% fileName='isEEInv_0_isEIInv_0_III_-1.0';
% fileName='isEEInv_1_isEIInv_1_pEE_1.0';
% fileName='isEEInv_1_isEIInv_1_pEI_1.0';
% fileName='isEEInv_1_isEIInv_1';
% fileName='isEEInv_0_isEIInv_0_ENonVar_1.0';
% fileName='Gamma_isEEInv_1_isEIInv_0_MeanI_90';
% fileName='Gamma_isEEInv_0_isEIInv_0_MeanI_270';
% fileName='Gamma_FixedESpk_isEEInv_1_isEIInv_1_pEI_0.0';
% fileName='Gamma_isEEInv_0_isEIInv_0_n_E_10';
% fileName='Gamma_tol_spk_0.01_isEEInvSuddenAP_1_isEIInvSuddenAP_0_NoUnrConn';
% fileName='Gamma_satisfycond_tol_spk_0.01';
% fileName='Non_Gamma_satisfycond_tol_spk_0.01';
% fileName='Gamma_satisfycond_tol_spk_0.0';
% fileName='Non_Gamma_satisfycond_tol_spk_0.0';
% fileName='pureInh';
% fileName='Non_Gamma_pureInh';
% fileName='Gamma_notsatisfycond';
% fileName='Non_Gamma_notsatisfycond';

fileName='tmp';

%% Which ones shall we record?
isRecord = 0;
isShowPatternOrig = 1;
isSavePatternOrig = 0;
isShowPattern = 0;
isShowTopology = 0;
isShowExcArrival = 0;
isShowInhArrival = 0;
isSaveFigure = 0;

%% The period of pattern that is normalized by the gamma period (30 ms).
Gamma_period=1/53*1e+3; % [ms]

%% To select whether we have the suddenly spike mechanism
isEEInvSuddenAP=0;pEEInvSuddenAP=0;
isEIInvSuddenAP=0;pEIInvSuddenAP=0;  % If isXXInvSuddenAP==1, pXXInvSuddenAP is the probability that the inputs invoke spike. pXXInvSuddenAP==1.0, all arrival inputs hit at the spike time of the postsynaptic neurons. pXXInvSuddenAP==0.0, no input hits at the spike time of the postsynaptic neurons.

%% The number of neurons in the network.
N=1000; % Osc. is ordered from 1 to 100.
% N=2;

%% Excitatory index.
N_ex=201; % A neuron is excitatory iff its index is greater or equal N_ex.
% N_ex=2;

%% The inh. input should arrive at tInhMin<=phase<=tInhMax. The exc. input should arrive at tExcMin<=phase<=tExcMax 
% tEEMeanDelay=2.5;   tEEMinDelayRange=-2.0;  tEEMaxDelayRange=30.0;  isEERandomlySelected=1.0;% [ms]
tEEMeanDelay=2.5;   tEEMinDelayRange=-0.3;  tEEMaxDelayRange=0.3;  isEERandomlySelected=1.0;    % [ms]
% tEEMeanDelay=3;   tEEMinDelayRange=-0.3;  tEEMaxDelayRange=0.3;  isEERandomlySelected=1.0;    % [ms]

tEIMeanDelay=1.3;   tEIMinDelayRange=-0.3;  tEIMaxDelayRange=0.3;  % [ms]
% tEIMeanDelay=1.3;   tEIMinDelayRange=-0.3;  tEIMaxDelayRange=30;  % [ms]

tIEMeanDelay=0.95;  tIEMinDelayRange=-0.1;  tIEMaxDelayRange=0.1;  % [ms]
% tIEMeanDelay=0.95;  tIEMinDelayRange=-0.1;  tIEMaxDelayRange=30;  % [ms]

tIIMeanDelay=0.6;   tIIMinDelayRange=-0.1;  tIIMaxDelayRange=0.1;  % [ms]
% tIIMeanDelay=0.6;   tIIMinDelayRange=-0.1;  tIIMaxDelayRange=30;  % [ms]

%% Alignment factor of the spike times
% ENonVar=1.0;EAlignment=-8.712*ENonVar+8.8;IAlignment=1e-6;
% EAlignment=-8.712*pEIInvSuddenAP+8.8;IAlignment=1e-6;
% tmpexp=pEIInvSuddenAP*2;EAlignment=1/(10^tmpexp);IAlignment=1e-6;
% EAlignment=2.0*1e-0;IAlignment=1e-6;
% EAlignment=1.0*1e+0;IAlignment=5.0*1e0;
% EAlignment=1*1e+2;IAlignment=1e+10;
% EAlignment=2*1e-1;IAlignment=1e-6;
EAlignment=1e-6;IAlignment=1e-6;

MeanI=88.405175905003787;
% MeanI=200;
n_E=-1;
MeanE=88.405175905003787+n_E*tEIMeanDelay*360/Gamma_period;

%% The number of presynaptic neurons contacts with a postsynaptic neuron.
pEE=0.02;
pEI=0.1;
% pEE=0;
% pEI=0;

pIE=0.3;
pII=0.1;    % Actually, it's 0.008 but because we have a small number and we want to investigate the I->I connections. Then, we set a bit higher.
% pIE=1.0;
% pII=1.0;    % Actually, it's 0.008 but because we have a small number and we want to investigate the I->I connections. Then, we set a bit higher.


%% Refractoriness
Refractory_period=0.5; % [ms]. Usually the value is 1.0. But because the delay is so short and we want to test for very align spike times, the I->I connections are not possible.

%% Parameters of LIF oscillator.
% gamma=1;I=1.6;
gamma=1;I=2.0;
% gamma=1;I=1.2;

%% Period of each oscillator.
T=1;

%% Upper and lower bound of weights.
lbEE=0.0001;ubEE=Inf;
lbEI=0.0001;ubEI=Inf;

lbIE=-Inf;ubIE=-0.0001;
% lbIE=-Inf;ubIE=-0.001;

lbII=-Inf;ubII=-0.0001;
% lbII=-Inf;ubII=-1;

%% Prohibit accidentally spiking before we want.
tol_boundary=0.05;

%% To make sure that the excitatory inputs bring the state more than one. This option is only used when the suddenly spking mechanism.
tol_spk=0.01;

%% The Inh. and Exc. inputs are not allowed to come within []
tSpkBoundary=1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataMode=2;
if (tSpkBoundary>=Refractory_period/Gamma_period)
    display('tSpkBoundary should be less than Refractory_period/Gamma_period');
    finish();
end

display('[Generate spiking times]');
spk_times=genSpikeTime(N,N_ex,EAlignment,IAlignment,MeanE,MeanI); % [degree] from 0 to 360.

%% The number of presynaptic neurons contacts with a postsynaptic neuron.
nEE=ceil(pEE*(N-N_ex+1));
nEI=ceil(pEI*(N-N_ex+1));
nIE=ceil(pIE*(N_ex-1));
nII=ceil(pII*(N_ex-1));

display(strcat('(E)->E:',num2str(nEE)));
display(strcat('(I)->E:',num2str(nIE)));
display(strcat('(E)->I:',num2str(nEI)));    
display(strcat('(I)->I:',num2str(nII)));       


if ((nEE<1) || (nEI<1) || (nIE<1) || (nII<1))
    display('#Connections is less than 1');
end

%% The inh. input should arrive at tInhMin<=phase<=tInhMax. The exc. input should arrive at tExcMin<=phase<=tExcMax 
tEEMin=Refractory_period/Gamma_period;tEEMax=1.0;
tEIMin=Refractory_period/Gamma_period;tEIMax=1.0;
tIEMin=Refractory_period/Gamma_period;tIEMax=1.0;
tIIMin=Refractory_period/Gamma_period;tIIMax=1.0;

tEEMinDelay=(tEEMeanDelay+tEEMinDelayRange)/Gamma_period;tEEMaxDelay=(tEEMeanDelay+tEEMaxDelayRange)/Gamma_period;
tEIMinDelay=(tEIMeanDelay+tEIMinDelayRange)/Gamma_period;tEIMaxDelay=(tEIMeanDelay+tEIMaxDelayRange)/Gamma_period;
tIEMinDelay=(tIEMeanDelay+tIEMinDelayRange)/Gamma_period;tIEMaxDelay=(tIEMeanDelay+tIEMaxDelayRange)/Gamma_period;
tIIMinDelay=(tIIMeanDelay+tIIMinDelayRange)/Gamma_period;tIIMaxDelay=(tIIMeanDelay+tIIMaxDelayRange)/Gamma_period;

osc_types=ones(1,N);osc_types(1,1:1:(N_ex-1))=-1; % Type of oscillator: -1 (inhibitory) and 1 (excitatory).

Natural_period=(-1/gamma)*log((1-I/gamma)/(-I/gamma));

periods=ones(1,N).*T;
gamma_i=ones(1,N).*gamma;
I_i=ones(1,N).*I;

spk_times=spk_times.*Gamma_period./360; % [ms]
statSpk_Times(spk_times,N_ex,N,Gamma_period);
spk_times=spk_times/Gamma_period; %[-] normalized by the gamma period
if (isShowPatternOrig == 1)
    % Show spike times
    figure(1)
    hold on
    
    for i=1:N
        if (osc_types(1,i)==-1)
            line([spk_times(:,i) spk_times(:,i)],[i i+0.5],'LineStyle','-','LineWidth',4,'Color',[0 0 1])
        else
            line([spk_times(:,i) spk_times(:,i)],[i i+0.5],'LineStyle','-','LineWidth',4,'Color',[1 0 0])
        end        
    end
    
    set(gca,'XTick',[0 1],'XTickLabel',{'';''});
    set(gca,'YTick',[0],'YTickLabel',{''});
    
    xlim([0 1]);
%     xlabel('t/T0');
%     grid on
    box on

    ylim([0 N+1]);
    
    if (isSavePatternOrig==1)
        print(strcat(fileName,'.eps'), '-depsc2');
    end
end

%% Generate connection matrix from osc. j (column) to osc. i (row).
display('[Determine whether two neurons are connected]');
conns=genConnMatrix(N,0);   % conns(i,j)==1 iff there is a connection from osc. j to i. conns(i,j)==-1 iff there is NO connection from osc. j to i.

%% Generate delay for postsynaptic neurons suddenly spike
display('[Generate delays]');

[delays,isESuddenAP,isISuddenAP]=genDelaysV4(conns,spk_times,N_ex,T,tol_eq, ...
                                            tEEMin, tEEMax, tEEMinDelay, tEEMaxDelay, ...
                                            tEIMin, tEIMax, tEIMinDelay, tEIMaxDelay, ...
                                            tIIMin, tIIMax, tIIMinDelay, tIIMaxDelay, ...
                                            tIEMin, tIEMax, tIEMinDelay, tIEMaxDelay, ...
                                            tSpkBoundary, ...
                                            isEEInvSuddenAP, isEIInvSuddenAP, ...
                                            pEEInvSuddenAP, pEIInvSuddenAP);    % delays(i,j)==NaN iff there is NO connection from osc. j to i.

% Get rid of unreasonable connections.
[conns,delays,isESuddenAP,isISuddenAP]=killConns(conns,delays,nEE,nEI,nIE,nII,N_ex,isESuddenAP,isISuddenAP,isEERandomlySelected);

if (isDelayCorrect(spk_times,periods,delays,osc_types,N,N_ex,tol_eq, ...
        tEEMin, tEEMax, tEEMinDelay, tEEMaxDelay, ...
        tEIMin, tEIMax, tEIMinDelay, tEIMaxDelay, ...
        tIIMin, tIIMax, tIIMinDelay, tIIMaxDelay, ...
        tIEMin, tIEMax, tIEMinDelay, tIEMaxDelay, ...
        conns, isESuddenAP, isISuddenAP, ...
        nEE, nEI, nIE, nII)==0)
    display('Incorrect delay');
    return;
end

% % Get delays that provide the same dynamics as reference
% delays=getDelayWithSameRelativeArr(N,spk_times,T,tol_eq,arr4eachOsc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isShowTopology == 1)
    showNetworkTopology(spk_times,T,conns,delays,N_ex);
end

lbInh=conns; % From I->X
ubInh=conns; % From I->X
lbExc=conns; % From E->X
ubExc=conns; % From E->X

lbInh(N_ex:N,1:(N_ex-1))=lbInh(N_ex:N,1:(N_ex-1)).*lbIE;
lbInh(1:(N_ex-1),1:(N_ex-1))=lbInh(1:(N_ex-1),1:(N_ex-1)).*lbII;
ubInh(N_ex:N,1:(N_ex-1))=ubInh(N_ex:N,1:(N_ex-1)).*ubIE;
ubInh(1:(N_ex-1),1:(N_ex-1))=ubInh(1:(N_ex-1),1:(N_ex-1)).*ubII;

lbExc(N_ex:N,N_ex:N)=lbExc(N_ex:N,N_ex:N).*lbEE;
lbExc(1:(N_ex-1),N_ex:N)=lbExc(1:(N_ex-1),N_ex:N).*lbEI;
ubExc(N_ex:N,N_ex:N)=ubExc(N_ex:N,N_ex:N).*ubEE;
ubExc(1:(N_ex-1),N_ex:N)=ubExc(1:(N_ex-1),N_ex:N).*ubEI;

% Get weights.
display('[Calculating weight]');
[weights,arr4eachOsc,v0,ret_msg]=getWeightsV4(size(spk_times,2),spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary,tol_spk,isESuddenAP,isISuddenAP);

if (isWeightCorrect(weights,osc_types,N,tol_eq,ret_msg)==0)
    return
end

% Save calculated data
if (isRecord == 1)
    display('[Recording the data]');
    save(strcat('../Data/',fileName,'.mat'),'N','tol_spk','lbInh','ubInh','lbExc','ubExc','conns','tol_boundary', ...
        'Gamma_period','N_ex','fileName','Rnd_Number','NoConn', ...
        'neuronModel','weights','arr4eachOsc','delays','v0', ...
        'spk_times','tBegin','tEnd','dt','dataMode','T','periods', ...
        'gamma_i','I_i','tol_eq','osc_types', ...
        'isEEInvSuddenAP','isEIInvSuddenAP', ...
        'EAlignment','IAlignment', ...
        'MeanI','MeanE', ...
        'pEE','pEI','pIE','pII', ...
        'Refractory_period', ...
        'tEEMeanDelay','tEEMinDelayRange','tEEMaxDelayRange','tEIMeanDelay','tEIMinDelayRange','tEIMaxDelayRange','tIIMeanDelay','tIIMinDelayRange','tIIMaxDelayRange','tIEMeanDelay','tIEMinDelayRange','tIEMaxDelayRange', ...
        'pEEInvSuddenAP','pEIInvSuddenAP');
    mysave(eclipse_ws,fileName,N,NoConn,weights,delays,v0,spk_times,T,tol_eq,gamma_i,I_i,Refractory_period,Gamma_period);
end

display('[Statistical information of the network]');
display(strcat('T0:',num2str(Natural_period),'=',num2str(Natural_period.*Gamma_period),' ms'));
statConnections(N,N_ex,conns);
statDelays(delays,N_ex,Gamma_period,tEEMeanDelay,tEEMinDelayRange,tEEMaxDelayRange,tEIMeanDelay,tEIMinDelayRange,tEIMaxDelayRange,tIIMeanDelay,tIIMinDelayRange,tIIMaxDelayRange,tIEMeanDelay,tIEMinDelayRange,tIEMaxDelayRange);
statWeights(weights,N_ex,Gamma_period);

if (strcmp(ret_msg,'Succes!!!')==1)    
    if (isShowPattern == 1)
        display('[Showing the pattern]');
        
        figure(3)
        dV=rand(1,N)*0; % perturb the initial voltages.
        tol_eq=1e-14;
        
        simNetworkEB(weights, delays, v0, spk_times, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, tol_eq, osc_types, dV, isSaveFigure, fileName, isShowExcArrival, isShowInhArrival);
%         maximize('all');
    end
end

end

function mysave(eclipse_ws,fileName,N,NoConn,weights,delays,v0,spk_times,T,tol_eq,gamma_i,I_i,Refractory_period,Gamma_period)
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
fwrite(fid, Refractory_period, 'double');
fwrite(fid, Gamma_period, 'double');

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
                %                 if (cmp(r,0,tol_eq)==0)
                if (cmp(r,0,2.220446049250314e-16)==0)
                    r=T;
                end
                return;
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

display(strcat('E: avg=',num2str(mean(spk_times(:,N_ex:N))),' ms, min=',num2str(min(spk_times(:,N_ex:N))),' ms, max=',num2str(max(spk_times(:,N_ex:N)))));
display(strcat('I: avg=',num2str(mean(spk_times(:,1:(N_ex-1)))),' ms, min=',num2str(min(spk_times(:,1:(N_ex-1)))),' ms, max=',num2str(max(spk_times(:,1:(N_ex-1))))));

end

function statWeights(weights,N_ex,Gamma_period)
E_E=0;E_E_n=0;
E_I=0;E_I_n=0;
I_E=0;I_E_n=0;
I_I=0;I_I_n=0;

N=size(weights,2);

for i=1:1:N
    for j=1:1:N
        if ~isnan(weights(i,j))
            if ((i>=N_ex) && (j>=N_ex))
                % E to E
                E_E_n=E_E_n+1;
                E_E(E_E_n)=weights(i,j);
            elseif ((i>=N_ex) && (j<N_ex))
                % I to E
                I_E_n=I_E_n+1;
                I_E(I_E_n)=weights(i,j);
            elseif ((i<N_ex) && (j>=N_ex))
                % E to I
                E_I_n=E_I_n+1;
                E_I(E_I_n)=weights(i,j);
            elseif ((i<N_ex) && (j<N_ex))
                % I to I
                I_I_n=I_I_n+1;
                I_I(I_I_n)=weights(i,j);
            end
        end
    end
end

% mean(E_E)
% mean(E_I)

display(strcat('E->E:',num2str(mean(E_E)),' std:',num2str(std(E_E))));
display(strcat('I->E:',num2str(mean(I_E)),' std:',num2str(std(I_E))));
display(strcat('E->I:',num2str(mean(E_I)),' std:',num2str(std(E_I))));
display(strcat('I->I:',num2str(mean(I_I)),' std:',num2str(std(I_I))));
end

function statDelays (delays,N_ex,Gamma_period, ...
                        tEEMeanDelay,tEEMinDelayRange,tEEMaxDelayRange, ...
                        tEIMeanDelay,tEIMinDelayRange,tEIMaxDelayRange, ...
                        tIIMeanDelay,tIIMinDelayRange,tIIMaxDelayRange, ...
                        tIEMeanDelay,tIEMinDelayRange,tIEMaxDelayRange)

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
% 
% mean(E_E)
% mean(E_I)

display(strcat('E->E:',num2str(mean(E_E)*Gamma_period),'/',num2str(mean(tEEMeanDelay)),' ms std:',num2str(std(E_E)*Gamma_period),' ms'));
display(strcat('I->E:',num2str(mean(I_E)*Gamma_period),'/',num2str(mean(tIEMeanDelay)),' ms std:',num2str(std(I_E)*Gamma_period),' ms'));
display(strcat('E->I:',num2str(mean(E_I)*Gamma_period),'/',num2str(mean(tEIMeanDelay)),' ms std:',num2str(std(E_I)*Gamma_period),' ms'));
display(strcat('I->I:',num2str(mean(I_I)*Gamma_period),'/',num2str(mean(tIIMeanDelay)),' ms std:',num2str(std(I_I)*Gamma_period),' ms'));
end

function val=isInOpenInterval(x_begin,x,x_end,tol_eq)
val=((cmp(x_begin,x,tol_eq)<0) && (cmp(x,x_end,tol_eq)<0));
end

function val=isInCloseInterval(x_begin,x,x_end,tol_eq)
val=((cmp(x_begin,x,tol_eq)<=0) && (cmp(x,x_end,tol_eq)<=0));
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

function showNetworkTopology(spk_times,T,conns,delays,N_ex)
noOscs=size(spk_times,2);

figure(10)
% Plot E->E
subplot(2,2,1);hold on
plot(spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');

for post=N_ex:noOscs    
    for pre=N_ex:noOscs
        if (conns(post,pre)==1)
            arriv=spk_times(1,pre)+delays(post,pre);
            plot([spk_times(1,pre) arriv],[pre post],'r');            
        end
    end
end

plot(T+spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');
title('E->E');
xlim([0 2*T]);

% Plot E->I
subplot(2,2,2);hold on
plot(spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');
plot(spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');

for post=1:(N_ex-1)
    for pre=N_ex:noOscs
        if (conns(post,pre)==1)
            arriv=spk_times(1,pre)+delays(post,pre);
            plot([spk_times(1,pre) arriv],[pre post],'r');            
        end
    end
end

plot(T+spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');
plot(T+spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');
title('E->I');
xlim([0 T/2]);

% Plot I->E
subplot(2,2,3);hold on
plot(spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');
plot(spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');

for post=N_ex:noOscs
    for pre=1:(N_ex-1)
        if (conns(post,pre)==1)
            arriv=spk_times(1,pre)+delays(post,pre);
            plot([spk_times(1,pre) arriv],[pre post],'b');            
        end
    end
end

plot(T+spk_times(1,N_ex:noOscs),N_ex:noOscs,'r*');
plot(T+spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');
title('I->E');
xlim([0 T/2]);

% Plot I->I
subplot(2,2,4);hold on
plot(spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');

for post=1:(N_ex-1)
    for pre=1:(N_ex-1)
        if (conns(post,pre)==1)
            arriv=spk_times(1,pre)+delays(post,pre);
            plot([spk_times(1,pre) arriv],[pre post],'b');            
        end
    end
end

plot(T+spk_times(1,1:(N_ex-1)),1:(N_ex-1),'b*');
title('I->I');
xlim([0 T/2]);

end

function statConnections(N,N_ex,conns)
    % The number of E->E outputs (1->M)
    % The number of E->I outputs (1->M)
    % The number of I->E outputs (1->M)
    % The number of I->I outputs (1->M)
    EE_output=zeros(N,1);    
    EI_output=zeros(N,1);
    IE_output=zeros(N,1);    
    II_output=zeros(N,1);

    % The number of E->E inputs (N->1)
    % The number of E->I inputs (N->1)
    % The number of I->E inputs (N->1)
    % The number of I->I inputs (N->1)    
    EE_input=zeros(N,1);    
    EI_input=zeros(N,1);
    IE_input=zeros(N,1);    
    II_input=zeros(N,1);    
    
    for j=1:1:N
        for i=1:1:N
            if (conns(i,j)==1)
                if ((N_ex<=j) && (N_ex<=i))
                    % E->E
                    EE_output(j,1)=EE_output(j,1)+1;
                    EE_input(i,1)=EE_input(i,1)+1;
                elseif ((N_ex<=j) && (i<=(N_ex-1)))
                    % E->I
                    EI_output(j,1)=EI_output(j,1)+1;
                    EI_input(i,1)=EI_input(i,1)+1;                    
                elseif ((j<=(N_ex-1)) && (N_ex<=i))
                    % I->E
                    IE_output(j,1)=IE_output(j,1)+1;
                    IE_input(i,1)=IE_input(i,1)+1;                    
                else
                    % I->I
                    II_output(j,1)=II_output(j,1)+1;
                    II_input(i,1)=II_input(i,1)+1;                    
                end
            end
        end
    end
    
    min_EE_output=min(EE_output(N_ex:N,1));
    max_EE_output=max(EE_output(N_ex:N,1));
    avg_EE_output=sum(EE_output(N_ex:N,1))/(N-N_ex+1);
    
    min_EI_output=min(EI_output(N_ex:N,1));
    max_EI_output=max(EI_output(N_ex:N,1));
    avg_EI_output=sum(EI_output(N_ex:N,1))/(N-N_ex+1);
    
    min_IE_output=min(IE_output(1:(N_ex-1),1));
    max_IE_output=max(IE_output(1:(N_ex-1),1));
    avg_IE_output=sum(IE_output(1:(N_ex-1):N,1))/(N_ex-1);

    min_II_output=min(II_output(1:(N_ex-1),1));
    max_II_output=max(II_output(1:(N_ex-1),1));
    avg_II_output=sum(II_output(1:(N_ex-1):N,1))/(N_ex-1);
    
    min_EE_input=min(EE_input(N_ex:N,1));
    max_EE_input=max(EE_input(N_ex:N,1));
    avg_EE_input=sum(EE_input(N_ex:N,1))/(N-N_ex+1);

    min_EI_input=min(EI_input(1:(N_ex-1),1));
    max_EI_input=max(EI_input(1:(N_ex-1),1));
    avg_EI_input=sum(EI_input(1:(N_ex-1),1))/(N_ex-1);    
    
    min_IE_input=min(IE_input(N_ex:N,1));
    max_IE_input=max(IE_input(N_ex:N,1));
    avg_IE_input=sum(IE_input(N_ex:N,1))/(N-N_ex+1);

    min_II_input=min(II_input(1:(N_ex-1),1));
    max_II_input=max(II_input(1:(N_ex-1),1));
    avg_II_input=sum(II_input(1:(N_ex-1),1))/(N_ex-1);        
    
%     display(strcat('E->N (E):',num2str(min_EE_output),'/',num2str(avg_EE_output),'/',num2str(max_EE_output)));
%     display(strcat('E->N (I):',num2str(min_EI_output),'/',num2str(avg_EI_output),'/',num2str(max_EI_output)));
%     display(strcat('I->N (E):',num2str(min_IE_output),'/',num2str(avg_IE_output),'/',num2str(max_IE_output)));    
%     display(strcat('I->N (I):',num2str(min_II_output),'/',num2str(avg_II_output),'/',num2str(max_II_output)));    
    
    display(strcat('(E)->E:',num2str(min_EE_input),'/',num2str(avg_EE_input),'/',num2str(max_EE_input)));
    display(strcat('(I)->E:',num2str(min_IE_input),'/',num2str(avg_IE_input),'/',num2str(max_IE_input)));
    display(strcat('(E)->I:',num2str(min_EI_input),'/',num2str(avg_EI_input),'/',num2str(max_EI_input)));    
    display(strcat('(I)->I:',num2str(min_II_input),'/',num2str(avg_II_input),'/',num2str(max_II_input)));       
end



function val=isDelayCorrect(spk_times,periods,delays,osc_types,N,N_ex,tol_eq, ...
    tEEMin, tEEMax, tEEMinDelay, tEEMaxDelay, ...
    tEIMin, tEIMax, tEIMinDelay, tEIMaxDelay, ...
    tIIMin, tIIMax, tIIMinDelay, tIIMaxDelay, ...
    tIEMin, tIEMax, tIEMinDelay, tIEMaxDelay, ...
    conns, isESuddenAP, isISuddenAP, ...
    nEE, nEI, nIE, nII)

% % E->E
% for i=N_ex:N    
%     if (sum(sum(~isnan(delays(i,N_ex:N))))~=nEE)
%         val=0;
% %         display(strcat('Osc. ',num2str(i),' E->E is not correct'));
% %         return;
%     end
% end
% 
% % I->E
% for i=N_ex:N    
%     if (sum(sum(~isnan(delays(i,1:(N_ex-1)))))~=nIE)
%         val=0;
% %         display(strcat('Osc. ',num2str(i),' I->E is not correct'));
% %         return;
%     end
% end
% 
% % E->I
% for i=1:(N_ex-1)
%     if (sum(sum(~isnan(delays(i,N_ex:N))))~=nEI)
%         val=0;
% %         display(strcat('Osc. ',num2str(i),' E->I is not correct'));
% %         return;
%     end
% end
% 
% % I->I
% for i=1:(N_ex-1)
%     if (sum(sum(~isnan(delays(i,1:(N_ex-1)))))~=nII)
%         val=0;
% %         display(strcat('Osc. ',num2str(i),' I->I is not correct'));
% %         return;
%     end
% end

for j=1:1:N
    for i=1:1:N
        if (conns(i,j)==1)
            if ((osc_types(1,j)==1) && (osc_types(1,i)==1))
                tMin=tEEMin;
                tMinDelay=tEEMinDelay;
                tMax=tEEMax;
                tMaxDelay=tEEMaxDelay;
            elseif ((osc_types(1,j)==1) && (osc_types(1,i)==-1))
                tMin=tEIMin;
                tMinDelay=tEIMinDelay;
                tMax=tEIMax;
                tMaxDelay=tEIMaxDelay;
            elseif ((osc_types(1,j)==-1) && (osc_types(1,i)==-1))
                tMin=tIIMin;
                tMinDelay=tIIMinDelay;
                tMax=tIIMax;
                tMaxDelay=tIIMaxDelay;
            elseif ((osc_types(1,j)==-1) && (osc_types(1,i)==1))
                tMin=tIEMin;
                tMinDelay=tIEMinDelay;
                tMax=tIEMax;
                tMaxDelay=tIEMaxDelay;
            end
            
            if ~((tMinDelay<=delays(i,j)) && (delays(i,j)<=tMaxDelay))
                val=0;
            else
                val=1;
            end
            
%             if      (((1==osc_types(j)) && ( 1==osc_types(i)) && (1==isESuddenAP(i,j))) || ...
%                     (( 1==osc_types(j)) && (-1==osc_types(i)) && (1==isISuddenAP(i,j))))
%                 ab_arrival_time=spk_times(1,j)+delays(i,j);
%                 [n,r]=mMod(spk_times(1,i),periods(1,i),ab_arrival_time,tol_eq);
% 
%                 if ((cmp(r,0.0,1e-15) == 0 ) || (cmp(r,1.0,1e-15) == 0 ))
%                     val=1;
%                 else
%                     val=0;
%                 end
%             else
%                 if ((osc_types(1,j)==1) && (osc_types(1,i)==1))
%                     tMin=tEEMin;
%                     tMinDelay=tEEMinDelay;
%                     tMax=tEEMax;
%                     tMaxDelay=tEEMaxDelay;
%                 elseif ((osc_types(1,j)==1) && (osc_types(1,i)==-1))
%                     tMin=tEIMin;
%                     tMinDelay=tEIMinDelay;
%                     tMax=tEIMax;
%                     tMaxDelay=tEIMaxDelay;
%                 elseif ((osc_types(1,j)==-1) && (osc_types(1,i)==-1))
%                     tMin=tIIMin;
%                     tMinDelay=tIIMinDelay;
%                     tMax=tIIMax;
%                     tMaxDelay=tIIMaxDelay;
%                 elseif ((osc_types(1,j)==-1) && (osc_types(1,i)==1))
%                     tMin=tIEMin;
%                     tMinDelay=tIEMinDelay;
%                     tMax=tIEMax;
%                     tMaxDelay=tIEMaxDelay;
%                 end
% 
%                 ab_arrival_time=spk_times(1,j)+delays(i,j);
%                 [n,r]=mMod(spk_times(1,i),periods(1,i),ab_arrival_time,tol_eq);
%                 val=(isInCloseInterval(tMinDelay,delays(i,j),tMaxDelay,tol_eq) && isInCloseInterval(tMin,r,tMax,tol_eq));                            
%             end                
            
            if (val==0)
                display(strcat('Delay from Osc.',num2str(j),' to Osc.',num2str(i),' is not correct'));
                return;
            end
        end
    end
end
end