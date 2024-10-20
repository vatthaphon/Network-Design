% clc;
clear all;close all

addpath('../Conn_matrix');
addpath('../Check_calculated_weights');
addpath('../Extract_pdf_from_the_paper');
addpath('../Gen_delays');
addpath('../Get_weights');
addpath('../LyapExp');
addpath('../Sim_network');
addpath('../Network_analysis');
addpath('../Data');

format long

names={ ...
%     'Gamma_params_tol_boundary_0_2_tol_spk_0_0_IncSpkVoltage_0_01' ...
%     'Gamma_params_tol_boundary_0_2_tol_spk_0_0_IncSpkVoltage_0_02' ...
%     'Gamma_params_tol_boundary_0_2_tol_spk_0_0_IncSpkVoltage_0_03' ...
%     'Gamma_params_tol_boundary_0_2_tol_spk_0_0_IncSpkVoltage_0_04' ...
%     'Gamma_params_tol_boundary_0_2_tol_spk_0_0_IncSpkVoltage_0_05' ...
%     'Gamma_notsatisfycond_params' ...
%     'Gamma_params_onlyInh' ...
    'Gamma_satisfycond_params' ...
%     'testE_E_NonSudden' ...
%     'testE_E_Sudden' ...
%     'testI_E_Non_Sudden' ...
%     'testI_E_Sudden' ...
%     'testI_I' ...
    };

tol_eq=eps;
for name=names
    filename=strcat(char(name),'.mat');    
    load(filename);
    
%     seed=3;
%     factor=1e-10;    
%     iter=1;
%     iterBegin=1;
% 
%     rand('seed',seed);    
%     init_perturbation=(rand(1,N)-rand(1,N))*factor;
%     [lambda,alpha]=maxLyapExpReturnMapV4(weights,delays,v0,iter,iterBegin,spk_times,N,T,tol_eq,gamma_i,I_i,init_perturbation,NoConn);

    
    rand('seed',101);
    init_perturbation=(rand(1,N)-rand(1,N))*1e-10;
    [lambda,alpha]=maxLyapExpReturnMapV4(weights,delays,v0,1,1,spk_times,N,T,tol_eq,gamma_i,I_i,init_perturbation,NoConn);            
%     plot(lambda)
%     log(1/init_perturbation(1,2))

%     nEvents=1;
%     noOscs=size(weights,1);
%     rand('seed',2);
%     init_perturbation=(rand(1,N)-rand(1,N))*1e-10;
%     v0;
%     v0=v0+init_perturbation;
%     pert_v0=nextSpkTimeV5(weights,delays,v0,spk_times,noOscs,tol_eq,gamma_i,I_i,NoConn,nEvents,T);
%     
%     sum(sum(abs(pert_v0-v0)))
end    

