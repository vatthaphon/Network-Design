function delays=genDelaysV2(conns,spk_times,N_ex,T,tol_eq)
N=size(spk_times,2);
delays=NaN(N);

% Set delays for excitatory connections
for j=N_ex:1:N
    for i=1:1:N
        if (conns(i,j)==1)            
            delays(i,j)=spk_times(1,i)-spk_times(1,j);
            if (cmp(delays(i,j),0,tol_eq)<=0)
                delays(i,j)=delays(i,j)+T;
                
                x=rand;
                while ((cmp(x, 0.9, tol_eq)>0) || (cmp(x, 0.1, tol_eq)<0))
                    x=rand;
                end
                
                delays(i,j)=delays(i,j)*x;
            else
                x=rand;
                while ((cmp(x, 0.9, tol_eq)>0) || (cmp(x, 0.1, tol_eq)<0))
                    x=rand;
                end
                
                delays(i,j)=delays(i,j)*(x/2+0.5);                
            end
            
        end
    end
end

% Set delays for inhibitory connections
for j=1:1:(N_ex-1)
    for i=1:1:N
        if (conns(i,j)==1)
            delays(i,j)=spk_times(1,i)-spk_times(1,j);
            if (cmp(delays(i,j),0,tol_eq)<=0)
                delays(i,j)=delays(i,j)+T;

                x=rand;
                while ((cmp(x, 0.9, tol_eq)>0) || (cmp(x, 0.1, tol_eq)<0))
                    x=rand;
                end

                delays(i,j)=delays(i,j)*x/2;                
            else
                x=rand;
                while ((cmp(x, 0.9, tol_eq)>0) || (cmp(x, 0.1, tol_eq)<0))
                    x=rand;
                end
                
                delays(i,j)=delays(i,j)+T*x/2;
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