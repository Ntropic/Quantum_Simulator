%second_test_CreateControlledCompleteSet_Cn_U2.m 
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

n=3;

[PU_gates,Cn]=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,0,0],1);

for i=1:length(PU_gates(:))
    fprintf(['\n<strong>Gate:</strong>  ' PU_gates(i).names ' \n']);
    fock_index=PU_gates(i).fock_index;
    fprintf(['<strong>Indexes:</strong>  [' num2str(fock_index(1)) ', ' num2str(fock_index(2)) '] \n']);
    FockPrint(PU_gates(i).matrix);
end