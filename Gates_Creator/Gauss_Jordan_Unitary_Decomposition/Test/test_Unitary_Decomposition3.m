%test_Unitary_Decomposition3.m
%Test Unitary Decomposition and Compare with shortened/optimized
%decomposition
clear all;
close all;
clc;

do_print=1;
do_tex=0;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Number of qubits
n=3;
sizle=2^n;

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

%Create/Load Controlled Unitaries
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);

%Create randum unitary matrix 
matrix=RandU(sizle)

%Transform
mode='opt';            %Optimized
[gates_opt,check_matrix]=Unitary_Decomposition(matrix,PU_gates,mode,'U_n','U_n',do_print);
mode='def';            %Non optimized - Default
[gates_def,check_matrix]=Unitary_Decomposition(matrix,PU_gates,mode,'U_n','U_n',do_print);

if do_tex==1
    %Create TeX File from gates variable
    file_name=[gates_opt.names '_optimized_short_15'];
    mode='expansion';
    depth=1;
    sep_dist=50;
    file_names=Gates2Tex(gates_opt,elem_gates,Comp_Gate_Merger(comp_gates,PU_gates),file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end

    %Create TeX File from gates variable
    file_name=[gates_def.names '_default'];
    mode='expansion';
    depth=2;
    sep_dist=50;
    file_names=Gates2Tex(gates_def,elem_gates,Comp_Gate_Merger(comp_gates,PU_gates),file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end
end