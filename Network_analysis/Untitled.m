t=linspace(0,100,100);
tau_d=3;
tau_r=0.1;
val=(1/(tau_d-tau_r))*(exp(-t/tau_d)-exp(-t/tau_r));

plot(t,val);

