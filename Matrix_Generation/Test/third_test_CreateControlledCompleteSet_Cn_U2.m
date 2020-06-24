%third_test_CreateControlledCompleteSet_Cn_U2.m 
%Create tex and pdf and image of the gates
clc;
clear all;
close all;

do_matrix=1;
do_tex=0;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

n=3;
m=5;
[PU_gates,Cn]=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0],1);

gates_def=PU_gates(m);


%Check matrix
if do_matrix==1
    [ matrix anc_matrix gates_expanded gates_indexes_expanded] = Gate2Matrix( elem_gates, comp_gates, gates_def, 2 )
end

%Create TeX File from gates variable
if do_tex==1  
    file_name=[gates_def.names '_expanded'];
    mode='expansion';
    depth=1;
    sep_dist=80;
    file_names=Gates2Tex(gates_def,elem_gates,Comp_Gate_Merger(comp_gates,PU_gates),file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end
end

