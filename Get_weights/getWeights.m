function [weights,arr4eachOsc,v0,ret_msg] = getWeights(spk_times, gamma_i, I_i, delays, T, ex_neuron_i, tol_eq, lbInh, ubInh, lbExc, ubExc, tol_boundary)
% GETWEIGHTS Determine the weight of the connections
%   [weights,ret_flag] = getWeights(spk_times, gamma_i, I_i, delays, T, ex_neuron_i, tol_eq)
% Input
%   spk_times   [ms] a row vector sorted from low to high: spk_times(1,1)=0;
%   gamma_i     a row vector containing decay rate of oscillator i.
%   I_i         a row vector containing constant current injected to oscillator i.
%   delays      [ms] a delay time matrix. The delay can be 0. When there is
%               no connection from osc. j (column) to osc. i (row), value of
%               the delay is NaN.
%   T           [ms] period of the pattern.
%   ex_neuron_i a number indicating that if the index of the oscillator is
%               greater than or equal 'ex_neuron_i', the oscillator is the
%               excitatory neuron.
%   tol_eq      two numbers are equal iff |x1-x2|<tol_eq.
% Output
%   weights     a weight matrix. w_ij represents a weight from osc. j (column) to
%               osc. i (row).
%   arr4eachOsc arrival times of inputs to each osc. with respect to the
%               spike time of each osc.
%   v0          a row vector, element of which is the voltage after spiking
%               of osc. 1.
%   ret_msg     returned message.
%   tol_boundary    to substract from the boundary of the condition to make
%                   sure that the inequality holds because Matlab uses the 
%                   equal sign in Ax<=b.
%
% Oscillators are indexed from 1 to size(spk_times,2). When the index of
% the oscillator is less than 'ex_neuron_i', the oscillator is an
% inhibitory neuron: otherwise, it is an excitatory neuron. Two inputs can
% arrive at a postsynaptic neuron at the same time!.
%
% Atthaphon Viriyopase, Donders Institute, 2011

nOscs = size(spk_times, 2);
weights=NaN(nOscs);
ret_msg='Succes!!!';

%% Check validity of ordering of spike times of the oscillators.
if (isConstraintSatisfy(spk_times,tol_eq) == 0)
    ret_msg='Spike times are not ordered from low to high';
    return
end

%% Determine absolute arrival of spikes
absolute_spk_arr_times=NaN(nOscs);
for i=1:1:nOscs
    absolute_spk_arr_times(i,:)=delays(i,:)+spk_times(1,:);
end

%% Check if a spike causes a postsynaptic neuron suddenly spike, i.e., the
%  relative arrival time equals spk_time.
relative_spk_arr_times=rem(absolute_spk_arr_times,T);
for i=1:1:nOscs
    for j=1:1:nOscs
        if ~isnan(relative_spk_arr_times(i,j))
            if (cmp(relative_spk_arr_times(i,j),spk_times(1,i),tol_eq) == 0)
                ret_msg='An input arrival time and spike time of the postsynaptic neuron are the same';
                return
            end
        end
    end
end

%% Determine arrival times of the input with respect to the spike time of
%  the postsynaptic neuron
arr_times_for_each_oscs=NaN(nOscs);
for i=1:1:nOscs
    for j=1:1:nOscs
        if ~isnan(absolute_spk_arr_times(i,j))
            if (cmp(spk_times(1,i),absolute_spk_arr_times(i,j),tol_eq)<0)
                arr_times_for_each_oscs(i,j)=rem(absolute_spk_arr_times(i,j)-spk_times(1,i),T);
            else
                tmp=spk_times(1,i)-absolute_spk_arr_times(i,j);
                arr_times_for_each_oscs(i,j)=T-tmp;
            end
        end        
    end
end

arr4eachOsc=arr_times_for_each_oscs;

%% Determine the weights between the oscillators
weights=NaN(nOscs);
for i=1:1:nOscs
    % Get the oscillators sending spikes to oscillator i
    from_oscs=find(~isnan(arr_times_for_each_oscs(i,:)));
    
    % Get the sorted arrival times of osc. i
    spk_arr_times=arr_times_for_each_oscs(i,~isnan(arr_times_for_each_oscs(i,:)));
    [spk_arr_times,Ix]=sort(spk_arr_times,'ascend');
    from_oscs=from_oscs(Ix);
    
    intv_btw_spks=calIntv_btw_Spks(spk_arr_times);
    
    % Index of the last spike arrival
    N=size(spk_arr_times,2)-1;
    
    % Form matrices for the optimization problem
    H=getH(N+1);
    F=getF(N+1);
    
    Aeq=exp(gamma_i(1,i).*spk_arr_times);
    beq=exp(gamma_i(1,i)*T)*(1-I_i(1,i)/gamma_i(1,i))+I_i(1,i)/gamma_i(1,i);    
    
    A=getA(N,gamma_i(1,i),intv_btw_spks);
    b=getb(N,gamma_i(1,i),intv_btw_spks,I_i(1,i),tol_boundary);
    
    A_dash=getA_dash(N,gamma_i(1,i),spk_arr_times);
    b_dash=getB_dash(N,gamma_i(1,i),I_i(1,i),spk_arr_times,tol_boundary);
    
    A=[A;A_dash];
    b=[b;b_dash];
    
    lb=getlb(from_oscs,ex_neuron_i,lbInh,lbExc);
    ub=getub(from_oscs,ex_neuron_i,ubInh,ubExc);
    
    
    options=optimset('display','off','LargeScale','off');
    
    [weight,fval,exitflag,output,lambda] = quadprog(H,F,A,b,Aeq,beq,lb,ub,[],options);
    display(strcat(printMsg(exitflag),': osc. ',num2str(i)));
    if (exitflag<=0)
        ret_msg=strcat(printMsg(exitflag),': osc. ',num2str(i));
    end
    
    weights(i,from_oscs)=weight';
end

%% Determine voltages of each oscillator after the first osc. spikes.
v0=NaN(1,nOscs);
v0(1,1)=0;
for i=2:1:nOscs
    duration=T-spk_times(1,i);
    [tSpan,V]=checkCalculatedWeights(arr4eachOsc(i,:),weights(i,:),duration,gamma_i(1,i),I_i(1,i));
    v0(1,i)=V(1,end);
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

function val=getub(from_oscs,ex_neuron_i, ubInh, ubExc)
n=size(from_oscs,2);

val=NaN(n,1);

for i=1:1:n
    if (from_oscs(1,i)>=ex_neuron_i)
        % Excitation
        val(i,1)=ubExc;
    else
        % Inhibition
        val(i,1)=ubInh;
    end    
end

end

function val=getlb(from_oscs,ex_neuron_i,lbInh,lbExc)
n=size(from_oscs,2);

val=NaN(n,1);

for i=1:1:n
    if (from_oscs(1,i)>=ex_neuron_i)
        % Excitation
        val(i,1)=lbExc;
    else
        % Inhibition
        val(i,1)=lbInh;
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
end

function val=getB_dash(N,gamma,I,spk_arr_times,tol_boundary)
val=zeros(N,1);
for j=1:1:N
    val(j,1)=(1-tol_boundary-I/gamma)*exp(gamma*spk_arr_times(1,j+1))+I/gamma;
end
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

function ret_flag=isConstraintSatisfy(spk_times,tol_eq)
ret_flag=1;

n=size(spk_times,2);

%% No two spike times that have the same value
for i=1:1:(n-1)
    if (cmp(spk_times(1,i), spk_times(1,i+1), tol_eq) > 0)
        ret_flag=0;
        return        
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