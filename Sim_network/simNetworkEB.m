function simNetworkEB(epsiMat, tauMat, initVoltages, tSpikes, tBegin, tEnd, dt, dataMode, T, gamma_i, I_i, tol_eq, osc_types, dV, isSaveFigure, fileName, isShowExcArrival, isShowInhArrival)
% epsiMat       strengths of the instantaneous synapse of the directed
%               graph. When there is no connection from osc. j (column) to osc. i (row),
%               tauMat(i,j)=NaN.
% tauMat        time delays. Its unit is the same as the unit of time. When
%               there is no connection from osc. j (column) to osc. i (row),
%               tauMat(i,j)=NaN.
% initVoltages  initial voltages of the oscillator. It's a row vector.
% tSpikes       List of spike pattern. It's a row vector. Its unit is the same
%               as the unit of time.
% tBegin        beginning of the simulation.
% tEnd          end of the simulation.
% dt            time step.
% dataMode      0=Raoul's data,1=Atthaphon's data
% T             period of the pattern. Its unit is the same as the unit of
%               time.
% gamma_i       a row vector containing gamma for each LIF.
% I_i           a row vector containing constant current.
% osc_types     a row vector containing -1=inh. neuron,1=exc. neuron
% dV            perturbation of initial voltage

global neuronModel
global NoConn

% 0: the Mirollo Strogatz neuron
% 1: the Leaky Integrate-and-Fire neuron
neuronModel=1;

tauMat(isnan(tauMat))=NoConn;

noOscs=size(initVoltages,2);

spkArrivalArrays=zeros(noOscs,1);
spkCSArrays=zeros(noOscs,1);
spkFromArrays=zeros(noOscs,1);
spkCounts=zeros(noOscs,1);

voltages=initVoltages+dV;
t=tBegin;
spikeTimes_i=0;
inhArrivalTime=zeros(noOscs,1);
inhArrivalTime_i=zeros(noOscs,1);
exArrivalTime=zeros(noOscs,1);
exArrivalTime_i=zeros(noOscs,1);

[results,spkCounts]=initialSpkArrivalArrays(spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,epsiMat,tauMat,tSpikes,dataMode,T);
spkArrivalArrays    =results(1:noOscs,:);
spkCSArrays         =results(noOscs+1:2*noOscs,:);
spkFromArrays       =results(2*noOscs+1:3*noOscs,:);

ref_Osc=1;
rel_V_i=0;                

while (cmp(t,tEnd,tol_eq)<=0)
    %% Check if there is a spike coming
    isSpikeComing=chkSpikeComing(t,noOscs,spkArrivalArrays,spkCounts,tol_eq);
    
    %% Get the synaptic weights
    if (isSpikeComing==1)
        [results,spkCounts,totSynStrengths,inhArrivalTime,inhArrivalTime_i,exArrivalTime,exArrivalTime_i]= ...
            getSynStrengths(t,noOscs,spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,inhArrivalTime,inhArrivalTime_i,exArrivalTime,exArrivalTime_i,tol_eq);        
        spkArrivalArrays    =results(1:noOscs,:);
        spkCSArrays         =results(noOscs+1:2*noOscs,:);
        spkFromArrays       =results(2*noOscs+1:3*noOscs,:);
        
%To be deleted
pre_volt=voltages;
        voltages=voltages+totSynStrengths';
    end
    
    %% Check if the neurons are spiking
    isSpiking=doesNeuronSpiking(voltages,tol_eq);
    
    %% Send spikes
    if sum(isSpiking>0)
        [results,spkCounts]=addSpkArrivalArrays(t,spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,epsiMat,tauMat,isSpiking);
        spkArrivalArrays    =results(1:noOscs,:);
        spkCSArrays         =results(noOscs+1:2*noOscs,:);
        spkFromArrays       =results(2*noOscs+1:3*noOscs,:);
    end
    
    %% Reset the voltages when they are spiking
    voltages=resetNeuronStates(voltages,isSpiking);
        
    %% Record spike times
    if sum(isSpiking>0)
        isSpikingNew=double(isSpiking);
        isSpikingNew(isSpiking==0)=NaN;
        spikeTimes_i=spikeTimes_i+1;
        spikeTimes(spikeTimes_i,:)=isSpikingNew.*t;   
        
        if (isSpiking(1,ref_Osc)==1)
            rel_V_i=rel_V_i+1;            
            rel_V(rel_V_i,:)=voltages;            
        end         
    end 
    
    %% Update the natural spiking time of oscillator
    oscsNaturalSpikeTime=calNatualSpikeTime(t,voltages,gamma_i,I_i);
    
    %% Update time by warping    
    nextEvent=getNextEvent(spkArrivalArrays,spkCounts,oscsNaturalSpikeTime,tEnd,tol_eq);
    interval=nextEvent-t;
    
    %% Update voltages because of time warping
    voltages=warpUpdate(voltages,interval,gamma_i,I_i);    
    t=nextEvent;
end

%% Raster plot
% subplot(2,1,1);
hold on
noOscs=size(initVoltages,2);
inhArrivalTime(inhArrivalTime==0)=NaN;
exArrivalTime(exArrivalTime==0)=NaN;
for i=1:1:noOscs
    if (osc_types(1,i)==-1)
        plot(spikeTimes(:,i),i,'b*','MarkerSize',12);
%         line([spikeTimes(:,i) spikeTimes(:,i)],[i i+0.5],'LineStyle','-','LineWidth',4,'Color',[0 0 1])
    else
        plot(spikeTimes(:,i),i,'r*','MarkerSize',12);
%         line([spikeTimes(:,i) spikeTimes(:,i)],[i i+0.5],'LineStyle','-','LineWidth',4,'Color',[1 0 0])
    end
        
    if (isShowExcArrival==1)
        plot(exArrivalTime(i,:),i,'ro','MarkerSize',2);        
    end
    
    if (isShowInhArrival==1)
        plot(inhArrivalTime(i,:),i,'b+','MarkerSize',2);        
    end    
end
box on;
% set(gca,'FontSize', 40);
set(gca,'XTick',[0 1],'XTickLabel',{'';''});
set(gca,'YTick',[0],'YTickLabel',{''});
ylim([0 noOscs+1]);
xlim([0 tEnd]);
% xlim([0.1 0.4]);
% xlabel('Time/T_0');
% ylabel('Osc. ID');
% title('Blue color=Inhibitory oscillator & Red color=Excitatory oscillator');

line([1 1],[0 noOscs+1],'LineStyle','-','LineWidth',0.1,'Color',[0 0 0])

% subplot(2,1,2);
% hold on
% % for i=1:1:noOscs
% 
% ex_osc_i=min(find(osc_types==1));
% inh_osc_i=max(find(osc_types==-1));
% 
% plot(rel_V(:,ex_osc_i),'-*r');
% plot(rel_V(:,inh_osc_i),'-*b');
% plot(rel_V(:,1),'-*k');
% % end
% xlabel('Time');
% ylabel('Relative voltage');

if (isSaveFigure==1)        
    print(strcat('../Data/',fileName,'.eps'),'-depsc2');
end

end

%% Initialize spike arrival arrays
function [val,spkCounts]=initialSpkArrivalArrays(spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,epsiMat,tauMat,tSpikes,dataMode,T)
global NoConn
if (dataMode==0)
    noOscs=size(epsiMat,1);
    notSpks=size(tSpikes,2);
    tmp=tauMat;
    tmp(tmp==-1)=NaN;
    tau=min(min(tmp));
    
    for i=2:1:notSpks
        if (cmp(tSpikes(1,i),tau,tol_eq)<0)
            isSpiking=zeros(1,noOscs);
            isSpiking(1,i)=1;
            
            [results,spkCounts]=addSpkArrivalArrays(tSpikes(1,i),spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,epsiMat,tauMat,isSpiking);
            spkArrivalArrays    =results(1:noOscs,:);
            spkCSArrays         =results(noOscs+1:2*noOscs,:);
            spkFromArrays       =results(2*noOscs+1:3*noOscs,:);
        else
            val=[spkArrivalArrays;spkCSArrays;spkFromArrays];
            return
        end
    end
elseif (dataMode==1)
    nOscs = size(tSpikes, 2);
    
    %% Add ghost inputs that cause by osc. pre
    for pre=1:1:nOscs
        for post=1:1:nOscs
            if (tauMat(post,pre)~=NoConn)
                if (pre==1)
                    t=tSpikes(1,pre);    % Because we assume that the simulation starts by the time that osc. 1 just spikes.
                else
                    t=tSpikes(1,pre)-T;
                end
                
                % Because the calculated initial voltage is just after the
                % event at t=0, we cannot use >=0.                
                while (cmp(t+tauMat(post,pre),0,eps)>0)
                    spkCounts(post,1)=spkCounts(post,1)+1;
                    
                    spkArrivalArrays(post,spkCounts(post,1))=t+tauMat(post,pre);
                    spkCSArrays(post,spkCounts(post,1))=epsiMat(post,pre);
                    spkFromArrays(post,spkCounts(post,1))=pre;
                    
                    % Sort the arrival time
                    for j=1:1:nOscs
                        [tmp,ix]=sort(spkArrivalArrays(j,1:spkCounts(j,1)));
                        spkArrivalArrays(j,1:spkCounts(j,1))=spkArrivalArrays(j,ix);
                        spkCSArrays(j,1:spkCounts(j,1))=spkCSArrays(j,ix);
                        spkFromArrays(j,1:spkCounts(j,1))=spkFromArrays(j,ix);
                    end
                    
                    t=t-T;
                end
                
            end
        end
    end
    
    val=[spkArrivalArrays;spkCSArrays;spkFromArrays];
elseif (dataMode==2)
    nOscs = size(tSpikes, 2);
    
    %% Add ghost inputs that cause by osc. pre
    for pre=1:1:nOscs
        for post=1:1:nOscs
            if (tauMat(post,pre)~=NoConn)
                t=tSpikes(1,pre)-T;
                
                % Because the calculated initial voltage is just after the
                % event at t=0, we cannot use >=0.
                while (cmp(t+tauMat(post,pre),0,eps)>0)
                    spkCounts(post,1)=spkCounts(post,1)+1;
                    
                    spkArrivalArrays(post,spkCounts(post,1))=t+tauMat(post,pre);
                    spkCSArrays(post,spkCounts(post,1))=epsiMat(post,pre);
                    spkFromArrays(post,spkCounts(post,1))=pre;
                    
                    % Sort the arrival time
                    for j=1:1:nOscs
                        [tmp,ix]=sort(spkArrivalArrays(j,1:spkCounts(j,1)));
                        spkArrivalArrays(j,1:spkCounts(j,1))=spkArrivalArrays(j,ix);
                        spkCSArrays(j,1:spkCounts(j,1))=spkCSArrays(j,ix);
                        spkFromArrays(j,1:spkCounts(j,1))=spkFromArrays(j,ix);
                    end
                    
                    t=t-T;
                end
                
            end
        end
    end
    
    val=[spkArrivalArrays;spkCSArrays;spkFromArrays];    
end
end

%% Check if there is a spike coming
function val=chkSpikeComing(t,noOscs,spkArrivalArrays,spkCounts,tol_eq)
    val=0;
    for i=1:1:noOscs
        if (spkCounts(i,1)>=1)
            if (cmp(spkArrivalArrays(i,1),t,tol_eq)<=0)
                val=1;
                return
            end
        end
    end    
end

%% Get the synaptic weights
function [results,spkCounts,totSynStrengths,inhArrivalTime,inhArrivalTime_i,exArrivalTime,exArrivalTime_i]=...
    getSynStrengths(t,noOscs,spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,inhArrivalTime,inhArrivalTime_i,exArrivalTime,exArrivalTime_i,tol_eq)
    
    totSynStrengths=zeros(noOscs,1);
    for i=1:1:noOscs
        if (spkCounts(i,1)>=1)
            if (cmp(spkArrivalArrays(i,1),t,tol_eq)<=0)
                spkArrival  =spkArrivalArrays(i,1:1:spkCounts(i,1));
                spkCS       =spkCSArrays(i,1:1:spkCounts(i,1));
                spkFrom     =spkFromArrays(i,1:1:spkCounts(i,1));
                                
                actI=(cmpMat(spkArrival,t,tol_eq)<=0);
                
                syns=spkCS(actI);                
                totSynStrengths(i,1)=sum(sum(syns));
                
                noActI=size(syns,2);
                                
                spkArrival(actI)=0;
                spkCS(actI)=0;
                spkFrom(actI)=0;                
                
                nonActI=xor(actI,1);
                
                spkArrivalArrays(i,1:1:spkCounts(i,1))  =[spkArrival(nonActI) spkArrival(actI)];
                spkCSArrays(i,1:1:spkCounts(i,1))       =[spkCS(nonActI) spkCS(actI)];
                spkFromArrays(i,1:1:spkCounts(i,1))     =[spkFrom(nonActI) spkFrom(actI)];
                spkCounts(i,1)                          =spkCounts(i,1)-noActI;       
                
                
                isExcitatory=sum(sum(cmp(syns,0,tol_eq)>0));  % When the coming spikes are excitatory, isExcitatory>0
                isInhibitory=sum(sum(cmp(syns,0,tol_eq)<0));  % When the coming spikes are inhibitory, isInhibitory>0
                
                if (isExcitatory>0)
                    exArrivalTime_i(i,1)=exArrivalTime_i(i,1)+1;
                    exArrivalTime(i,exArrivalTime_i(i,1))=t;
                end
                
                if (isInhibitory>0)
                    inhArrivalTime_i(i,1)=inhArrivalTime_i(i,1)+1;
                    inhArrivalTime(i,inhArrivalTime_i(i,1))=t;
                end
            end
        end
    end     
    results=[spkArrivalArrays;spkCSArrays;spkFromArrays];
end


%% Send spikes
function [val,spkCounts]=addSpkArrivalArrays(t,spkArrivalArrays,spkCSArrays,spkFromArrays,spkCounts,epsiMat,tauMat,isSpiking)
global NoConn

    connectedNodes=(tauMat~=NoConn);
    [row,col] = find(connectedNodes);
    
    n=size(row,1);
    for i=1:1:n
        pre=col(i);post=row(i);
        if (isSpiking(1,pre)==1)
            spkCounts(post,1)=spkCounts(post,1)+1;
            
            spkArrivalArrays(post,spkCounts(post,1))=t+tauMat(post,pre);            
            spkCSArrays(post,spkCounts(post,1))=epsiMat(post,pre);
            spkFromArrays(post,spkCounts(post,1))=pre;
        end
    end
    
    % Sort the arrival time
    n=size(spkArrivalArrays,1);
    for i=1:1:n
        [tmp,ix]=sort(spkArrivalArrays(i,1:spkCounts(i,1)));
        spkArrivalArrays(i,1:spkCounts(i,1))=spkArrivalArrays(i,ix);
        spkCSArrays(i,1:spkCounts(i,1))=spkCSArrays(i,ix);
        spkFromArrays(i,1:spkCounts(i,1))=spkFromArrays(i,ix);        
    end
    val=[spkArrivalArrays;spkCSArrays;spkFromArrays];
end

%% Reset the voltages when they are spiking
function val=resetNeuronStates(voltages,isSpiking)
global neuronModel

    if (neuronModel==0)
        val=resetMSStates(voltages,isSpiking);
    elseif (neuronModel==1)
        val=resetLIFStates(voltages,isSpiking);
    end 
end

function val=resetMSStates(voltages,isSpiking)
    voltages(isSpiking)=0;
    val=voltages;
end

function val=resetLIFStates(voltages,isSpiking)
    voltages(isSpiking)=0;
    val=voltages;
end

%% Check if the neurons are spiking
function val=doesNeuronSpiking(voltages,tol_eq)
global neuronModel

    if (neuronModel==0)
        val=doesMSSpiking(voltages,tol_eq);
    elseif (neuronModel==1)
        val=doesLIFSpiking(voltages,tol_eq);
    end    
end

function val=doesMSSpiking(voltages,tol_eq)
    v_max=1.0;
    
    val=(cmpMat(voltages,v_max,tol_eq)>=0);
end

function val=doesLIFSpiking(voltages,tol_eq)
    v_max=1.0;

    val=(cmpMat(voltages,v_max,tol_eq)>=0);
end

%% Real time that oscillators will spike
function oscsNaturalSpikeTime=calNatualSpikeTime(t,voltages,gamma_i,I_i)
global neuronModel

    if (neuronModel==0)
    elseif (neuronModel==1)
        oscsNaturalSpikeTime=calNatualSpikeTimeLIF(t,voltages,gamma_i,I_i);
    end
end

function oscsNaturalSpikeTime=calNatualSpikeTimeLIF(t,voltages,gamma_i,I_i)
    oscsNaturalSpikeTime=t+(-1./gamma_i.*log((1-I_i./gamma_i)./(voltages-I_i./gamma_i)));
end

%% Update voltages of the neuron by warping time
function voltages=warpUpdate(voltages,interval,gamma_i,I_i)
global neuronModel

    if (neuronModel==0)
    elseif (neuronModel==1)
        voltages=warpUpdateLIF(voltages,interval,gamma_i,I_i);
    end
end

function voltages=warpUpdateLIF(voltages,interval,gamma_i,I_i)
    voltages=exp(-gamma_i.*interval).*voltages+I_i./gamma_i-I_i./gamma_i.*exp(-gamma_i.*interval);
end

%% Update voltages of the neuron
function val=update(voltages,dt,gamma_i,I_i)
    val=voltages+F(voltages,gamma_i,I_i).*dt;
end

function val=F(voltages,gamma_i,I_i)
global neuronModel

    if (neuronModel==0)
        val=FMS(voltages);
    elseif (neuronModel==1)
        val=FLIF(voltages,gamma_i,I_i);
    end
end

function val=FMS(voltages)
    b=3;
    
    val=(1/b)*(exp(b)-1)*exp(-b*voltages);
end

function val=FLIF(voltages,gamma_i,I_i)
    val=I_i-gamma_i.*voltages;
end

%% Misc.
function nextEvent=getNextEvent(spkArrivalArrays,spkCounts,oscsNaturalSpikeTime,tEnd,tol_eq)
    nextEvent=tEnd+1;
    n=size(spkArrivalArrays,1);
    for i=1:1:n
        if (spkCounts(i,1)>=1)
            if (cmp(spkArrivalArrays(i,1),nextEvent,tol_eq)<0)
                nextEvent=spkArrivalArrays(i,1);
            end
        end
    end
    
    for i=1:1:n
        if (cmp(oscsNaturalSpikeTime(1,i),nextEvent,tol_eq)<0)
            nextEvent=oscsNaturalSpikeTime(1,i);
        end
    end
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

function true_time=conv2TrueTimeLIF(dimensionless_time)
    I=1.2;
    T=log(I/(I-1));
    
    true_time=dimensionless_time*T;
end

function val=cmpMat(x,y,tol_eq)
[n,m]=size(x);
val=NaN(n,m);

[o,p]=size(y);

for i=1:1:n
    for j=1:1:m
        if ((o==1) && (p==1))
            val(i,j)=cmp(x(i,j),y,tol_eq);
        else
            val(i,j)=cmp(x(i,j),y(i,j),tol_eq);
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