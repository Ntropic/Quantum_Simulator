%test_Expansion2Matrix.m
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));
%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')


%Prepare gates
gates_table=Gates_Tables_Prepare(comp_gates,elem_gates)

%Choose gate
i=12;

%Error model
t_step=0.001;
T1=1;
T2=1;

phi=0;
E=0.003;    %Error rate
delta=sqrt(10/3*E); 
E1=5/4*E;

errors=[t_step,T1,T2,phi,delta,E1];

%Create Innitial density matrix
wave=[0 1 0 0]'; % fock2vec(1,2,[1 0;0 1],[1 1i]/sqrt(2));
rho0=full(sparse(Wave2Density(wave,1)));
DensityBar3( rho0,'Coloring',true,'AbsValue',true);
[rho sizes anc_sizes] = Noisy_Expansion2Matrix(rho0,errors,gates_table(i));
DensityBar3( rho,'Coloring',true,'AbsValue',true);

%Calculate ideal density
m=gates_table(i).matrix;
rho_i=m*rho0*m';
%Calculate fidelity %But with a different distance measure 
error_prob=Hilbert_Schmidt(rho_i,rho)