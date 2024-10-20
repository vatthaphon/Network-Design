function [weights,arr4eachOsc,v0,ret_msg] = getWeightsV2(nOscs,spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary)
% GETWEIGHTS Determine the weight of the connections
%   [weights,ret_flag] = getWeightsV2(nOscs,spk_times,periods,gamma_i,I_i,osc_types,delays,tol_eq,lbInh,ubInh,lbExc,ubExc,tol_boundary)
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
% Atthaphon Viriyopase, Donders Institute, 2011

weights=NaN(nOscs);
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


%% Check if a spike causes a postsynaptic neuron suddenly spike, i.e., the
%  relative arrival time equals spk_time.
for i=1:1:nOscs
    for j=1:1:nOscs
        if ~isnan(relative_spk_arr_times(i,j))
            if (cmp(relative_spk_arr_times(i,j),periods(1,i),tol_eq) == 0)
                ret_msg='An input arrival time and spike time of the postsynaptic neuron are the same';
                return
            end
        end
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
    
    intv_btw_spks=calIntv_btw_Spks(spk_arr_times);
    
    % Index of the last spike arrival
    N=size(spk_arr_times,2)-1;
    
    % Form matrices for the optimization problem
    H=getH(N+1);
    F=getF(N+1);
    
    Aeq=exp(gamma_i(1,i).*spk_arr_times);
    beq=exp(gamma_i(1,i)*periods(1,i))*(1-I_i(1,i)/gamma_i(1,i))+I_i(1,i)/gamma_i(1,i);    
    
    A=getA(N,gamma_i(1,i),intv_btw_spks);
    b=getb(N,gamma_i(1,i),intv_btw_spks,I_i(1,i),tol_boundary);
    
    A_dash=getA_dash(N,gamma_i(1,i),spk_arr_times);
    b_dash=getB_dash(N,gamma_i(1,i),I_i(1,i),spk_arr_times,tol_boundary);

    A_lb_ub=getA_lb_ub(from_oscs,osc_types);
    b_lb_ub=getb_lb_ub(from_oscs,i,lbExc,ubInh,osc_types);    
    
    A=[A;A_dash;A_lb_ub];
    b=[b;b_dash;b_lb_ub];
        
%     options=optimset('display','off','LargeScale','off','Algorithm','interior-point-convex');
%     
%     [weight,fval,exitflag,output,lambda] = quadprog(H,F,A,b,Aeq,beq,[],[],[],options);
%     
%     if (exitflag<=0)
%         display(strcat('Try determineAnalytically mode for osc. ',num2str(i)));
%         [weight,ret_flag]=determineAnalytically(weight,spk_arr_times,from_oscs,i,osc_types,periods(1,i),gamma_i(1,i),I_i(1,i),ubInh,lbExc,tol_eq);
%         
%         if ((ret_flag==0) && (strcmp(ret_msg,'Succes!!!')==1))
%             ret_msg=strcat(printMsg(exitflag),': osc. ',num2str(i));
%         end
%     end    

    [weight,err,lm] = qpas(H,F,A,b,Aeq,beq,[],[],0);
    if (cmp(err,0,tol_eq)~=0)
        display(strcat('Osc. ',num2str(i),' do not obtain the optimal solution'));
        ret_msg='Not Succes!!!';        
        return
    end
    
    weights(i,from_oscs)=weight';
end

% %% Determine initial voltages of each oscillator
% v0=NaN(1,nOscs);
% for i=1:1:nOscs
%     duration=periods(1,i)-spk_times(1,i);
%     [tSpan,V]=checkCalculatedWeights(arr4eachOsc(i,:),weights(i,:),duration,gamma_i(1,i),I_i(1,i));
%     v0(1,i)=V(1,end);
% end

%% Determine initial voltages of each oscillator
v0=NaN(1,nOscs);
for i=1:1:nOscs
    duration=periods(1,i)-spk_times(1,i);
    if (i==3)
        i
    end
    v0(1,i)=findInitialV(arr4eachOsc(i,:),weights(i,:),duration,gamma_i(1,i),I_i(1,i),tol_eq);    
end

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

function [weight,ret_flag]=determineAnalytically(old_weight,spk_arr_times,from_oscs,to_osc,osc_types,period,gamma,I,ubInh,lbExc,tol_eq)
n=size(from_oscs,2);

weight=NaN(n,1);
for i=n:-1:2
    if (osc_types(1,from_oscs(1,i))==1)
        weight(i,1)=lbExc(to_osc,from_oscs(1,i));        
    else
        weight(i,1)=ubInh(to_osc,from_oscs(1,i));        
    end        
end

t=period;
v0=1;
for i=n:-1:2
    inv_t=spk_arr_times(1,i)-t;
    v0=warpUpdate(v0,inv_t,gamma,I);
    
    v0=v0-weight(i,1);
    
    if (cmp(v0,1,tol_eq)>=0)
        display('Upper bound of weight problem');
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
   display('Voltage reaches 1 before being influenced by an input');
    
   ret_flag=0;
   weight=old_weight;
   return    
end

weight(1,1)=v0-tmpv0;

if (((cmp(weight(1,1),ubInh(to_osc,from_oscs(1,1)),tol_eq)>0) && (osc_types(1,from_oscs(1,1))==-1)) || ...
    ((cmp(weight(1,1),lbExc(to_osc,from_oscs(1,1)),tol_eq)<0) && (osc_types(1,from_oscs(1,1))==1)))
   display(strcat('The first weight should be ',num2str(weight(1,1))));
    
   ret_flag=0;
   weight=old_weight;
   return
end

ret_flag=1;

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



function val=getb_lb_ub(from_oscs,to_osc,lbExc,ubInh,osc_types)
n=size(from_oscs,2);
val=zeros(n,1);
for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==1)
        % Exc. neuron
        val(i,1)=-lbExc(to_osc,from_oscs(1,i));
    else
        % Inh. neuron
        val(i,1)=ubInh(to_osc,from_oscs(1,i));
    end
end
end

function val=getA_lb_ub(from_oscs,osc_types)
n=size(from_oscs,2);
val=zeros(n);
for i=1:1:n
    if (osc_types(1,from_oscs(1,i))==1)
        % Exc. neuron
        val(i,i)=-1;
    else
        % Inh. neuron
        val(i,i)=1;
    end
end
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

function val=getSpk_arr_times(spk_arr_times,j)
% Input
%   j       range from 0 to N.
% Output
%   
val=spk_arr_times(1,j+1);
end

function val=getIntv_btw_spks(intv_btw_spks,p)
% GETINTV_BTW_SPKS Get intv_btw_spks at index p
%   val = getIntv_btw_spks(intv_btw_spks,p)
% Input
%   intv_btw_spks   a row vector containing duration between the consequtive 
%                   arrival times.
%   p               an index of the arrival times ranging from -1,0,...,N-1.
% Output
%   val             intv_btw_spk

val=intv_btw_spks(1,p+2);
end

function val=getb(N,gamma,intv_btw_spks,I,tol_boundary)
% GETB Get the b matrix for the optimization problem
%   val = getb(N,gamma,intv_btw_spks,I)
% Input
%   N           the last index the arrival times. Hence, the number of all
%               arrival times is N+1.
%   gamma       decay rate of current oscillator.
%   intv_btw_spks   a row vector containing duration between the consequtive 
%                   arrival times
%   I           the current to the current oscillator.
% Output
%   val         a b matrix, size of which is (N+1)x1.

b=NaN(N+1,1);

b(1,1)=(1-tol_boundary)-getZ(getIntv_btw_spks(intv_btw_spks,-1),I,gamma);
for j=1:1:N
    sum=0;
    for l=0:1:(j-1)
        sum=sum+getZ(getIntv_btw_spks(intv_btw_spks,l-1),I,gamma)*exp(-gamma*sumIntv_btw_Spks(intv_btw_spks,l,j-1));
    end
    b(j+1,1)=(1-tol_boundary)-getZ(getIntv_btw_spks(intv_btw_spks,j-1),I,gamma)-sum;
end
val=b;

end

function val=getZ(intv_btw_spk,I,gamma)
val=I/gamma-I/gamma*exp(-gamma*intv_btw_spk);
end

function val=getA_dash(N,gamma,spk_arr_times)
val=zeros(N,N+1);
for j=1:1:N
    for k=0:1:(j-1)
        val(j,k+1)=exp(gamma*spk_arr_times(1,k+1));
    end
end

val=[zeros(1,N+1);val];
end

function val=getB_dash(N,gamma,I,spk_arr_times,tol_boundary)
val=zeros(N,1);
for j=1:1:N
    val(j,1)=(1-tol_boundary-I/gamma)*exp(gamma*spk_arr_times(1,j+1))+I/gamma;
end

tmp=(1-tol_boundary)+I/gamma*exp(-gamma*spk_arr_times(1,1))-I/gamma;
val=[tmp;val];
end


function val=getA(N,gamma,intv_btw_spks)
% GETA Get the A matrix for the optimization problem
%   val = getA(N,gamma,intv_btw_spks)
% Input
%   N           the last index of the arrival times. Hence, the number of all
%               arrival times is N+1.
%   gamma       decay rate of current oscillator.
%   intv_btw_spks   a row vector containing duration between the consequtive 
%                   arrival times
% Output
%   val         an A matrix, size of which is (N+1)x(N+1).

val=eye(N+1);
    
for j=1:1:N
    for jj=0:1:(j-1)
        val(j+1,jj+1)=exp(-gamma*sumIntv_btw_Spks(intv_btw_spks,jj,j-1));
    end
end
    
end

function val=getF(n)
val=zeros(n,1);
end

function val=getH(n)
val=eye(n);
end

function val=sumIntv_btw_Spks(intv_btw_spks,i_b,i_e)
val=0;
for i=i_b:1:i_e
    val=val+getIntv_btw_spks(intv_btw_spks,i);
end

end

function val=calIntv_btw_Spks(spk_arr_times)
% CALINTV_BTW_SPKS Calculate intervals between spike arrival times
%   val = calIntv_btw_Spks(spk_arr_times)
% Input
%   spk_arr_times   a row vector containing sorted (from low to high) 
%                   arrival times of the inputs.
% Output
%   val             a row vector containing intervals between arrival of
%                   inputs. size(val,2)=size(spk_arr_times,2).
%                   val(1,1)=spk_arr_times(1,1);

n=size(spk_arr_times,2);
val=NaN(1,n);
N=n-1;

val(1,1)=getSpk_arr_times(spk_arr_times,0);
for p=0:1:(N-1)
    val(1,p+2)=getSpk_arr_times(spk_arr_times,p+1)-getSpk_arr_times(spk_arr_times,p);
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
                return
            end
            n=n+1;
        end
    end
end
end



