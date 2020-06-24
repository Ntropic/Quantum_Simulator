%pcolormat_test.m
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-2});
addpath(genpath(shortened));
%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

%Test
n=7;
s=ceil(log2(n+1));
H=Exchange_Hamiltonian(n);

mode={'fock_gray_tex',n,2};
PColorMat(H,mode,800,['Particle_Mode_N_2'],1260);