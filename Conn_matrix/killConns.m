function [conns,delays,isESuddenAP,isISuddenAP]=killConns(conns,delays,nEE,nEI,nIE,nII,N_ex,isESuddenAP,isISuddenAP,isEERandomlySelected)
N=size(conns,1);
for i=1:1:N
    for j=1:1:N
        if isnan(delays(i,j))
            conns(i,j)=0;
        end
    end
end

for i=1:1:N
    if (i<N_ex)
        %% From E->I
        true_n_EI=sum(sum(conns(i,N_ex:N)));
        if (true_n_EI<nEI)
            % If the number of the satisfied delays are less, we use all.
            display(strcat('The number of E->I condidates is too less:',num2str(i),':',num2str(true_n_EI),'/',num2str(nEI)));            
        else
            % Pick up nEI connections randomly.
            tmp=conns(i,N_ex:N).*rand(1,N-N_ex+1);
            [tmp,Ix]=sort(tmp,'descend');
            
            conns(i,(N_ex-1)+Ix(1,(nEI+1):end))=0;
            delays(i,(N_ex-1)+Ix(1,(nEI+1):end))=NaN;
            
            isISuddenAP(i,(N_ex-1)+Ix(1,(nEI+1):end))=NaN;
        end
               
        %% From I->I
        true_n_II=sum(sum(conns(i,1:(N_ex-1))));
        if (true_n_II<nII)
            display(strcat('The number of I->I condidates is too less:',num2str(i),':',num2str(true_n_II),'/',num2str(nII)));            
        else
            % Pick up nII connections randomly.
            tmp=conns(i,1:(N_ex-1)).*rand(1,N_ex-1);
            [tmp,Ix]=sort(tmp,'descend');
            
            conns(i,Ix(1,(nII+1):end))=0;
            delays(i,Ix(1,(nII+1):end))=NaN;
        end
        
    else
        %% From E->E
        true_n_EE=sum(sum(conns(i,N_ex:N)));
        if (true_n_EE<nEE)
            display(strcat('The number of E->E condidates is too less:',num2str(i),':',num2str(true_n_EE),'/',num2str(nEE)));            
        else
            if (isEERandomlySelected==1)
                % Pick up nEE connections randomly.
                tmp=conns(i,N_ex:N).*rand(1,N-N_ex+1);
                [tmp,Ix]=sort(tmp,'descend');
                
                conns(i,(N_ex-1)+Ix(1,(nEE+1):end))=0;
                delays(i,(N_ex-1)+Ix(1,(nEE+1):end))=NaN;
                
                isESuddenAP(i,(N_ex-1)+Ix(1,(nEE+1):end))=NaN;
            else
                % Pick up nEE the shortest connections.
                tmpDelay=delays(i,N_ex:N);
                [tmpDelay,Ix]=sort(tmpDelay,'ascend');
                
                conns(i,(N_ex-1)+Ix(1,(nEE+1):end))=0;
                delays(i,(N_ex-1)+Ix(1,(nEE+1):end))=NaN;
                
                isESuddenAP(i,(N_ex-1)+Ix(1,(nEE+1):end))=NaN;
            end
        end
                
        %% From I->E
        true_n_IE=sum(sum(conns(i,1:(N_ex-1))));
        if (true_n_IE<nIE)
            display(strcat('The number of I->E condidates is too less:',num2str(i),':',num2str(true_n_IE),'/',num2str(nIE)));
        else
            % Pick up nIE connections randomly.
            tmp=conns(i,1:(N_ex-1)).*rand(1,N_ex-1);
            [tmp,Ix]=sort(tmp,'descend');
            
            conns(i,Ix(1,(nIE+1):end))=0;
            delays(i,Ix(1,(nIE+1):end))=NaN;
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