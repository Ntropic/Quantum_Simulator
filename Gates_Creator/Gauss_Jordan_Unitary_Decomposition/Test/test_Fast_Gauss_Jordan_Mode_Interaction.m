%test_Fast_Gauss_Jordan_Interaction.m
clear all;
close all;
clc;

do_print=1;
do_tex=0;
do_anim=0;
do_gates=0;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

N=2;    %Photons
t=pi/4;

Hij=Gray_Exchange_Hamiltonian_Particles(N);
Uij=expm(-1i*Hij*t);
n=log2(size(Uij,1));

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
%Create/Load Controlled Unitaries for different sizes
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);
other_gates=Comp_Gate_Merger(comp_gates,PU_gates);

gates_opt=Fast_Gauss_Jordan_Decomposition(Uij,do_print,'U_n','U_n');

%Create unitary matrix 
if do_print==1
    FockPrint(Uij,'dense',[N,2])
end
if do_tex==1
    [filename]=Sparse2Tex(Hij,['MI_H_' num2str(N)],'tikz_gray_fock_dense_def_subs_var_width'); 
    tex2pdf2preview(filename,1,1000);
    [filename]=Sparse2Tex(Uij,['MI_U_' num2str(N) '_t_' num2str(t)],'unit_tikz_gray_fock_dense_def_subs_var_width'); 
    tex2pdf2preview(filename,1,1000);
end
if do_anim==1
    Unitary_Anim(Hij,[0,2*pi],{'fock',n,2,0,16},800,['Particle_Mode_N_' num2str(N)]);
end

%% Gates
if do_gates==1
    file_name=[gates_opt.names '_gate'];
    mode='expansion';
    depth=1;
    sep_dist=35;
    Gates2Circ2Preview(gates_opt,elem_gates,other_gates,['MI_F_GJ_' num2str(N)],mode,depth,sep_dist);
end