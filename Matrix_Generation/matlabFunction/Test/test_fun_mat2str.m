%test_fun_mat2str.m
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
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

test_gate=elem_gates(12).fun_mat

[ c_mat,c_vars ] = fun_mat2str( test_gate )
[ fun ] = str2fun_mat( c_mat,c_vars )

[ A_new ] = Embed_Gate( fun ,[1 3] ,4 )