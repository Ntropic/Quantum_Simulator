%test_Mode_Interaction_animation.m
% Animates the unitary evolution operator of the mode interaction in gray
% code numbering
clear all;
close all;
clc;

do_print=0;
do_tex=0;
do_anim=1;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Number of qubits per mode
N=4;
n=2*N;
sizle=2^n;

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

%Create/Load Controlled Unitaries for different sizes
H=Gray_Exchange_Hamiltonian(N);

%Create unitary matrix 
if do_print==1
    FockPrint(H,'dense',[2^N-1,2])
end
if do_tex==1
    [filename]=SparseFull2Tex(H,['Particle_Mode_N_' num2str(N)],'tikz_gray_fock_dense_def_subs_var_width'); 
    tex2pdf2preview(filename,1,1000);
end
if do_anim==1
    %mode={'fock_gray_mp4',2^N-1,2};
    mode={'none_gif',2^N-1,2};
    Unitary_Anim(H,[0,4*pi,800],mode,800,['Particle_Mode_N_' num2str(N)],1260);
end