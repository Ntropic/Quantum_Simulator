%test_Gauss_Jordan_RandU.m
clear all;
close all;
clc;

do_print=0;
do_tex=1;
do_anim=1;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

N=4;    %Qubits
t=pi/2;

%Hij=Gray_Exchange_Hamiltonian_Particles(N);
Uij=RandU(2^N);
n=log2(size(Uij,1));

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
%Create/Load Controlled Unitaries for different sizes
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);
other_gates=Comp_Gate_Merger(comp_gates,PU_gates);

[gates_opt,check_matrix]=Gauss_Jordan_Decomposition(Uij,PU_gates,'U_n','U_n',do_print);
max(check_matrix(:))

%Create unitary matrix 
% if do_tex==1
%     [filename]=Sparse2Tex(Hij,['Particle_Mode_N_' num2str(N)],'tikz_gray_fock_dense_def_subs_var_width'); 
%     tex2pdf2preview(filename,1,1000);
% end
% if do_anim==1
%      Unitary_Anim(Hij,[0,2*pi],{'fock',3,2,0,16},800,['Particle_Mode_N_' num2str(N)]);
% end

%% Gates

file_name=[gates_opt.names '_gate'];
mode='expansion';
depth=1;
sep_dist=25;
Gates2Circ2Preview(gates_opt,elem_gates,other_gates,['RandU_GJ_' num2str(n)],mode,depth,sep_dist);

