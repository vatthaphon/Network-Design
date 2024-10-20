function [weights,arr4eachOsc,v0,ret_msg] = getWeightsV3(nOscs,spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary)
% GETWEIGHTS Determine the weight of the connections
%   [weights,ret_flag] = getWeightsV3(nOscs,spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary)
% Input
%   nOscs       a number of oscillators.
%   spk_times   a row vector containing spike time of oscillators from 1
%               to nOscs.
%   periods     a row vector containing period of oscillators from 1 to
%               nOscs. !!! These periods are not natural period of osc. i.
%               !!!
%   gamma_i     a row vector containing decay rate of oscillator i.
%   I_i         a row vector containing constant current injected to oscillator i.
%   osc_types   a row vector containing either -1 (inhibitory) or 1 (excitatory).
%   delays      a delay time matrix. The delay can be 0. When there is
%               no connection from osc. j (column) to osc. i (row), value of
%               the delay is NaN.
%   tol_eq      two numbers are equal iff |x1-x2|<tol_eq.
%   lbInh           a lower bound matrix of inhibitory connection from osc.
%                   j (column) to osc. i (row). When there is no
%                   connection from osc. j to osc. i, value is NaN.
%   ubInh           an upper bound matrix of inhibitory connection from osc.
%                   j (column) to osc. i (row). Where thre is no connection
%                   , value is NaN.
%   lbExc           a lower bound matrix of excitatory connection from osc.
%                   j (column) to osc. i (row). When there is no
%                   connection from osc. j to osc. i, value is NaN.
%   ubExc           an upper bound matrix of excitatory connection from osc.
%                   j (column) to osc. i (row). Where thre is no connection
%                   , value is NaN.
%   tol_boundary    to make sure that the constraint is met. Substract from 
%                   the boundary of the condition to make sure that the 
%                   inequality holds because Matlab uses the equal sign in Ax<=b.
% Output
%   weights     a weight matrix. w_ij represents a weight from osc. j (column) to
%               osc. i (row).
%   arr4eachOsc arrival times of inputs to each osc. with respect to the
%               spike time of each osc.
%   v0          a row vector, element of which is the voltage after spiking
%               of osc. 1.
%   ret_msg     returned message.
%   tol_boundary    
% 
% Calculate the weights that make the postsynaptic neurons spike suddenly
% Atthaphon Viriyopase, Donders Institute, 2011

ret_msg='Succes!!!';

%% Determine absolute arrival of spikes.
absolute_spk_arr_times=NaN(nOscs);
for i=1:1:nOscs
    % All synaptic inputs arriving at osc. i.
    absolute_spk_arr_times(i,:)=delays(i,:)+spk_times(1,:);
end

%% Calculate arrival of spike times relative to the period of the
%  postsynaptic oscillator.
relative_spk_arr_times=NaN(nOscs);
for i=1:1:nOscs
    for j=1:1:nOscs
        [n,r]=mMod(spk_times(1,i),periods(1,i),absolute_spk_arr_times(i,j),tol_eq);
        relative_spk_arr_times(i,j)=r;
    end
end

arr4eachOsc=relative_spk_arr_times;

%% Determine the weights between the oscillators
weights=NaN(nOscs);
for i=1:1:nOscs
    % Get the oscillators sending spikes to oscillator i
    from_oscs=find(~isnan(relative_spk_arr_times(i,:)));
    
    % Get the sorted arrival times of osc. i
    spk_arr_times=relative_spk_arr_times(i,~isnan(relative_spk_arr_times(i,:)));
    [spk_arr_times,Ix]=sort(spk_arr_times,'ascend');
    from_oscs=from_oscs(Ix);
        
    % Index of the last spike arrival
    N=size(spk_arr_times,2)-1;
    
    % Form matrices for the optimization problem
    H=getH(N+1);
    F=getF(N+1);
    
    Aeq=exp(gamma_i(1,i).*spk_arr_times);
    beq=exp(gamma_i(1,i)*periods(1,i))*(1-I_i(1,i)/gamma_i(1,i))+I_i(1,i)/gamma_i(1,i);    
    
    A_dash=getA_dash(N,gamma_i(1,i),spk_arr_times,tol_eq);
    b_dash=getB_dash(N,gamma_i(1,i),I_i(1,i),spk_arr_times,tol_boundary);
    
    lb=getlb(from_oscs,i,lbInh,lbExc,osc_types);
    ub=getub(from_oscs,i,ubInh,ubExc,osc_types);
                
    options=optimset('display','off','LargeScale','off');    
    [weight,fval,exitflag,output,lambda] = quadprog(H,F,A_dash,b_dash,Aeq,beq,lb,ub,[],options);

    if (exitflag<=0)
        display(strcat('Try determineAnalytically mode for osc. ',num2str(i)));
        [weight,ret_flag]=determineAnalytically(weight,spk_arr_times,from_oscs,i,osc_types,periods(1,i),gamma_i(1,i),I_i(1,i),ubInh,tol_eq);
        
        if ((ret_flag==0) && (strcmp(ret_msg,'Succes!!!')==1))
            ret_msg=strcat(printMsg(exitflag),': osc. ',num2str(i));
        end
    end
    
    weights(i,from_oscs)=weight';
end

%% Determine initial voltages of each oscillator
v0=NaN(1,nOscs);
for i=1:1:nOscs
    duration=periods(1,i)-spk_times(1,i);

%     [tSpan,V]=checkCalculatedWeights(arr4eachOsc(i,:),weights(i,:),duration,gamma_i(1,i),I_i(1,i));
%     % Voltage of osc. just after the duration.
%     v0(1,i)=V(1,end);

    v0(1,i)=findInitialV(arr4eachOsc(i,:),weights(i,:),duration,gamma_i(1,i),I_i(1,i),tol_eq);    
end

end

function [weight,ret_flag]=determineAnalytically(old_weight,spk_arr_times,from_oscs,to_osc,osc_types,period,gamma,I,ubInh,tol_eq)
% Check if all inputs are inhibitory
n=size(from_oscs,2);
for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==1)
        display('I give up to calculate the weights');
        ret_flag=0;
        weight=old_weight;
        return
    end
end

weight=NaN(n,1);
for i=n:-1:2
    weight(i,1)=ubInh(to_osc,from_oscs(1,i));
end

t=period;
v0=1;
for i=n:-1:2
    inv_t=spk_arr_times(1,i)-t;
    v0=warpUpdate(v0,inv_t,gamma,I);
    
    v0=v0-weight(i,1);
    
    if (cmp(v0,1,tol_eq)>=0)
        display('To large upper bound of inh. neuron');
        ret_flag=0;
        weight=old_weight;
        return        
    end
    
    t=t-abs(inv_t);
end

inv_t=spk_arr_times(1,1)-t;
v0=warpUpdate(v0,inv_t,gamma,I);
tmpv0=warpUpdate(0,spk_arr_times(1,1),gamma,I);

if (cmp(tmpv0,1,tol_eq)>=0)
   ret_flag=0;
   weight=old_weight;
   display('Voltage reaches 1 before being influenced by an input');
   return    
end

weight(1,1)=v0-tmpv0;

if (cmp(weight(1,1),ubInh(to_osc,from_oscs(1,1)),tol_eq)>0)
   ret_flag=0;
   weight=old_weight;
   display(strcat('The first weight should be',num2str(weight(1,1))));
   return
end

ret_flag=1;

end

function v0=findInitialV(arr4eachOsc,weights,duration,gamma,I,tol_eq)
% arr4eachOsc   a row vector containing arrival times for each oscillator.
% weights       a row vector containing weights that arrive each time.
% duration      tEnd

arr4eachOsc=arr4eachOsc(~isnan(arr4eachOsc));
weights=weights(~isnan(weights));

[arr4eachOsc,ix]=sort(arr4eachOsc,'ascend');
weights=weights(ix);

% Get arrival times before duration.
n=size(arr4eachOsc,2);
valid_i=0;
for i=1:1:n
    if (cmp(arr4eachOsc(1,i),duration,tol_eq)<=0)
        validArr4eachOsc(1,valid_i+1)=arr4eachOsc(1,i);
        validWeights(1,valid_i+1)=weights(1,i);
        valid_i=valid_i+1;
    end    
end

% v0 at duration.
v0=0;
t=0;
if (valid_i>0)
    n=size(validArr4eachOsc,2);
    for i=1:1:n
        dt=validArr4eachOsc(1,i)-t;
        
        v0=warpUpdate(v0,dt,gamma,I);
        
        v0=v0+validWeights(1,i);
        
        t=t+dt;
    end
    v0=warpUpdate(v0,duration-validArr4eachOsc(1,n),gamma,I);
else
    v0=warpUpdate(v0,duration,gamma,I);
end



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

function val=printMsg(exitflag)
    if (exitflag==1)
        val='Function converged to a solution x';
    elseif (exitflag==3)
        val='Change in the objective function value was smaller than the specified tolerance';
    elseif (exitflag==4)
        val='Local minimizer was found';
    elseif (exitflag==0)        
        val='Number of iterations exceeded options.MaxIter';
    elseif (exitflag==-2)            
        val='Problem is infeasible';
    elseif (exitflag==-3)        
        val='Problem is unbounded';
    elseif (exitflag==-4)
        val='Current search direction was not a direction of descent. No further progress could be made';
    elseif (exitflag==-7)
        val='Magnitude of search direction became too small. No further progress could be made';
    end
end

function A_lb_ub=getA_lb_ub(from_oscs,osc_types)
n=size(from_oscs,2);
A_lb_ub=zeros(n);

for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==-1)
        % Inhibition
        A_lb_ub(i,i)=1;
    else
        % Excitation
        A_lb_ub(i,i)=-1;
    end            
end

end

function b_lb_ub=getb_lb_ub(from_oscs,to_osc,ubInh,lbExc,osc_types)
n=size(from_oscs,2);
b_lb_ub=zeros(n,1);

for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==-1)
        % Inhibition
        b_lb_ub(i,1)=ubInh(to_osc,from_oscs(1,i));
    else
        % Excitation
        b_lb_ub(i,1)=-lbExc(to_osc,from_oscs(1,i));
    end          
end
end

function val=getub(from_oscs,to_osc,ubInh,ubExc,osc_types)
n=size(from_oscs,2);

val=NaN(n,1);

for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==-1)
        % Inhibition
        val(i,1)=ubInh(to_osc,from_oscs(1,i));        
    else
        % Excitation
        val(i,1)=ubExc(to_osc,from_oscs(1,i));                
    end    
end

end

function val=getlb(from_oscs,to_osc,lbInh,lbExc,osc_types)
n=size(from_oscs,2);

val=NaN(n,1);

for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==-1)
        % Inhibition
        val(i,1)=lbInh(to_osc,from_oscs(1,i));        
    else
        % Excitation
        val(i,1)=lbExc(to_osc,from_oscs(1,i));                
    end
end

end

function val=getA_dash(N,gamma,spk_arr_times,tol_eq)
val=zeros(N+1,N+1);

for i=2:1:N+1
    t=spk_arr_times(1,i);
    
    j=1;
    while (cmp(spk_arr_times(1,j),t,tol_eq)<0)
        val(i,j)=exp(gamma*spk_arr_times(1,j));
        j=j+1;
    end
end

end

function val=getB_dash(N,gamma,I,spk_arr_times,tol_boundary)
val=zeros(N,1);
for j=1:1:N
    val(j,1)=(1-tol_boundary-I/gamma)*exp(gamma*spk_arr_times(1,j+1))+I/gamma;
end

tmp=(1-tol_boundary)+I/gamma*exp(-gamma*spk_arr_times(1,1))-I/gamma;
val=[tmp;val];
end

function val=getF(n)
val=zeros(n,1);
end

function val=getH(n)
val=eye(n);
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



