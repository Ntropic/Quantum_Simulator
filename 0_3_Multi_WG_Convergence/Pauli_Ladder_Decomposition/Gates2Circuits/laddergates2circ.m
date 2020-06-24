%laddergates2circ.m
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
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
load('../clifford_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,clifford_gates);




for i=[1,3,5,7,9]
    gates_opt=clifford_gates(i);
    file_name=[gates_opt.names 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=25;

    file_names=Gates2Tex(gates_opt,elem_gates,other_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
end

for i=[1,3,5]
    gates_opt=comp_gates(i);
    file_name=[gates_opt.names 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=25;

    file_names=Gates2Tex(gates_opt,elem_gates,other_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
end