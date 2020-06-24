%test_Unitary_Decomposition_Sparse.m
%Test Unitary Decomposition for Sparse Matrices
clear all;
close all;
clc;

do_print=1;
do_tex=0;
do_anim=0;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Number of qubits per mode
n=2;
s=ceil(log2(n+1));
t=pi/2;

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

%Create/Load Controlled Unitaries for different sizes
PU_gates=CreateCompleteSet_Cn_U2(2*s-1,elem_gates,comp_gates,[1,1,0]);

H=Gray_Exchange_Hamiltonian_Particles(n);
matrix=expm(-1i*H*t);
%matrix=RandU(4);

[gates_opt,check_matrix]=Gauss_Jordan_Decomposition(matrix,PU_gates,'U_n','U_n',do_print);

%Create unitary matrix 
if do_print==1
    FockPrint(H,'dense',[n,2])
end
if do_tex==1
    [filename]=Sparse2Tex(H,['Particle_Mode_N_' num2str(N)],'tikz_gray_fock_dense_def_subs_var_width'); 
    tex2pdf2preview(filename,1,1000);
end
if do_anim==1
    Unitary_Anim(H,[0,2*pi],{'gray_fock',n,s,0,16},800,['Particle_Mode_N_' num2str(n)]);
end