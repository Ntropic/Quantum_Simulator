%XY_quantum_gates.m
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
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

%-------------------------------------------------------------------------

name='XX_YY';
size=2;
anc_size=0;
gate_steps={'XX','YY'};
index_steps={[1 2],[1 2]};
param_steps={[],[]};
gate_string='e^{-i\frac{\phi}{2}(x x+y y)}';


gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
FockPrint(gate.matrix)

filename=Save_Gate_list('XY_quantum_gates','xy_gates',gate);
fprintf('File has been saved');
