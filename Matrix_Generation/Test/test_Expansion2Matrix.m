%test_Expansion2Matrix.m
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
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

i=12; %3;

%Prepare gates
gates_table=Gates_Tables_Prepare(comp_gates,elem_gates)

%Create Matrix for gate i

[gates_matrix sizes anc_sizes] = Expansion2Matrix(gates_table(i));

full(gates_matrix)