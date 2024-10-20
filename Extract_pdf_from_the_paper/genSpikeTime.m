function val=genSpikeTime(N,N_ex,width_PC,width_BC,MeanE,MeanI)
% pdf of excitatory neurons
mu_PC=MeanE; % degree in [0,360] of the gamma oscillation that E neuron spike times concentrate.
sigma_PC=8.812368160574666*width_PC;
% sigma_PC=width_PC;

% pdf of inhibitory neurons
mu_BC=MeanI; % degree in [0,360] of the gamma oscillation that I neuron spike times concentrate.
sigma_BC=0.727776975869855*width_BC;

val=NaN(1,N);

% Inhibition (Inhibitory cell)
for i=1:1:(N_ex-1)
    val(1,i)=-1;
    
    while ((val(1,i)<=0) || (val(1,i)>=360))
        val(1,i) = rem(abs(norminv(rand,mu_BC,sigma_BC)),360);
    end
end

% Excitation (Pyrimadal cell)
for i=N_ex:1:(N)
    val(1,i)=-1;
    
    while ((val(1,i)<=0) || (val(1,i)>=360))
        val(1,i) = rem(abs(norminv(rand,mu_PC,sigma_PC)),360);
    end    
end

end