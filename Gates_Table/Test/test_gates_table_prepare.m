%Circ_XY_All_Interaction.m
clear all;
close all;
clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));
%Load gates
load('../elem_gates.mat','-mat')
load('../comp_gates.mat','-mat')

gates_table=Gates_Tables_Prepare(comp_gates,elem_gates)

%Compare matrices
len=length(gates_table);

for i=1:len
    curr_name=gates_table(i).names;
    if isa(curr_name,'char')
        curr_name2=curr_name;
    else
        curr_name2=curr_name{1};
    end
    fprintf([curr_name2 '\n'])
    U_approx=full(Expansion2Matrix(gates_table(i)));
    
    if isa(U_approx,'sym')
        simplify(U_approx-gates_table(i).matrix,10)
    else
        U_approx-gates_table(i).matrix
    end
    
    fprintf('------------------------------------------\n')
end