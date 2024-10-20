function [delays,isESuddenAP,isISuddenAP]=genDelaysV4(conns, spk_times, N_ex, T, tol_eq, ...
                            tEEMin, tEEMax, tEEMinDelay, tEEMaxDelay, ...
                            tEIMin, tEIMax, tEIMinDelay, tEIMaxDelay, ...
                            tIIMin, tIIMax, tIIMinDelay, tIIMaxDelay, ...
                            tIEMin, tIEMax, tIEMinDelay, tIEMaxDelay, ...
                            tSpkBoundary, ...
                            isEEInvSuddenAP, isEIInvSuddenAP, ...
                            pEEInvSuddenAP, pEIInvSuddenAP)
N=size(spk_times,2);
delays=NaN(N);

% Set delays for inhibitory connections
for j=1:1:(N_ex-1)
    for i=1:1:N
        if (conns(i,j)==1)
            if (i>=N_ex)
                % I->E         
                ref_spk_time=spk_times(1,i)-T;
                delays(i,j)=calDelay(ref_spk_time,j,spk_times(1,j),i,spk_times(1,i),tIEMin,tIEMax,tIEMinDelay,tIEMaxDelay,tSpkBoundary,tol_eq,T);
            else
                % I->I
                ref_spk_time=spk_times(1,i)-T;
                delays(i,j)=calDelay(ref_spk_time,j,spk_times(1,j),i,spk_times(1,i),tIIMin,tIIMax,tIIMinDelay,tIIMaxDelay,tSpkBoundary,tol_eq,T);
            end            
        end
    end
end

% Set delays for excitatory connections
isESuddenAP=NaN(N,N);
isISuddenAP=NaN(N,N);
for j=N_ex:1:N
    for i=1:1:N
        if (conns(i,j)==1)   
            if (i>=N_ex)
                % E->E      
                if ((isEEInvSuddenAP==1) && (rand(1,1)<=pEEInvSuddenAP))
                    delays(i,j)=spk_times(1,i)-spk_times(1,j);
                    if (cmp(delays(i,j),0,tol_eq)<=0)
                        delays(i,j)=delays(i,j)+T;
                    end      
                    
                    if ~((tEEMinDelay<=delays(i,j)) && (delays(i,j)<=tEEMaxDelay))
                        delays(i,j)=NaN;
                        isESuddenAP(i,j)=NaN;
                    else
                        isESuddenAP(i,j)=1;                        
                    end
                else                    
                    ref_spk_time=spk_times(1,i)-T;
                    delays(i,j)=calDelay(ref_spk_time,j,spk_times(1,j),i,spk_times(1,i),tEEMin,tEEMax,tEEMinDelay,tEEMaxDelay,tSpkBoundary,tol_eq,T);
                    
                    if isnan(delays(i,j))
                        isESuddenAP(i,j)=NaN;
                    else
                        isESuddenAP(i,j)=0;
                    end
                end
            else
                % E->I
                if (isEIInvSuddenAP==1) && (rand(1,1)<=pEIInvSuddenAP)
                    delays(i,j)=spk_times(1,i)-spk_times(1,j);
                    if (cmp(delays(i,j),0,tol_eq)<=0)
                        delays(i,j)=delays(i,j)+T;
                    end      
                    
                    if ~((tEIMinDelay<=delays(i,j)) && (delays(i,j)<=tEIMaxDelay))
                        delays(i,j)=NaN;
                        isISuddenAP(i,j)=NaN;
                    else
                        isISuddenAP(i,j)=1;
                    end                    
                else
                    ref_spk_time=spk_times(1,i)-T;
                    delays(i,j)=calDelay(ref_spk_time,j,spk_times(1,j),i,spk_times(1,i),tEIMin,tEIMax,tEIMinDelay,tEIMaxDelay,tSpkBoundary,tol_eq,T);
                    
                    if isnan(delays(i,j))
                        isISuddenAP(i,j)=NaN;
                    else
                        isISuddenAP(i,j)=0;
                    end                    
                end
            end
        end
    end
end

end

function delay=calDelay(ref_spk_time,from,from_SpkTime,to,to_SpkTime,tXXMin,tXXMax,tXXMinDelay,tXXMaxDelay,tSpkBoundary,tol_eq,T)
tMinDelay=from_SpkTime+tXXMinDelay;
tMaxDelay=from_SpkTime+tXXMaxDelay;

cnt=0;
[tMin tMax]=getAllowableRange(ref_spk_time, T, tXXMin, tXXMax,tSpkBoundary,tol_eq);
while (isIntersect(tMin,tMax,tMinDelay,tMaxDelay,tol_eq)==0)
    if (cnt>10)
%         display(strcat(num2str(from),'->',num2str(to),' is not possible'));
        delay=NaN;
        return
    end
    
    ref_spk_time=ref_spk_time+T;
    [tMin tMax]=getAllowableRange(ref_spk_time, T, tXXMin, tXXMax,tSpkBoundary,tol_eq);
    
    cnt=cnt+1;
end

if (cmp(tMin,tMinDelay,tol_eq)<0)
    tMin=tMinDelay;
end

if (cmp(tMaxDelay,tMax,tol_eq)<0)
    tMax=tMaxDelay;
end

x=rand;
delay=(tMin-from_SpkTime)+(tMax-tMin)*x;
end

function val=isIntersect(tMin,tMax,tMinDelay,tMaxDelay,tol_eq)
    val=((isInOpenInterval(tMin,tMinDelay,tMax,tol_eq) || isInOpenInterval(tMin,tMaxDelay,tMax,tol_eq)) || ...
        (isInOpenInterval(tMinDelay,tMin,tMaxDelay,tol_eq) || isInOpenInterval(tMinDelay,tMax,tMaxDelay,tol_eq)));
end

function val=isInOpenInterval(x_begin,x,x_end,tol_eq)
    val=((cmp(x_begin,x,tol_eq)<0) && (cmp(x,x_end,tol_eq)<0));
end

function [tMin tMax]=getAllowableRange(ref_spk_time, T, tXXMin, tXXMax,tSpkBoundary,tol_eq)
if (cmp(ref_spk_time+tXXMin,ref_spk_time+tSpkBoundary,tol_eq)<=0)
    tMin=ref_spk_time+tSpkBoundary;
else
    tMin=ref_spk_time+tXXMin;
end

if (cmp(ref_spk_time+tXXMax,ref_spk_time+T-tSpkBoundary,tol_eq)<=0)
    tMax=ref_spk_time+tXXMax;
else
    tMax=ref_spk_time+T-tSpkBoundary;
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