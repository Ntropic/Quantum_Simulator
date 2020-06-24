%test_SparseFull2Tex.m
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

%Test CreateControlledCompleteSet_Cn_U2 for n
%Create/Load Controlled Unitaries
n=3;
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);
S=PU_gates(11).matrix;



[filename]=SparseFull2Tex(S,'Cn_test1','none_def_subs_var_width');
tex2pdf2preview(filename,1,600);
[filename]=SparseFull2Tex(S,'Cn_test2',{'tikz_fock_def_subs_var_width',1,n});
tex2pdf2preview(filename,1,600);
[filename]=SparseFull2Tex(S,'Cn_test3','tikz_fock_dense_def_subs_var_width'); %Best
tex2pdf2preview(filename,1,600);
[filename]=SparseFull2Tex(S,'Cn_test4','tikz_index_def_subs_var_width');
tex2pdf2preview(filename,1,600);