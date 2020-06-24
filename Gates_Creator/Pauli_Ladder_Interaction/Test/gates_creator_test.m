%gates_creator_test.m
clc;
clear all;
close all;

depth=3;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-depth});
addpath(genpath(shortened));

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
load('clifford_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

%-------------------------------------------------------------------------
n=4;
gate=Pauli_Ladder_Gates(n);
% gate2=Pauli_Ladder_Gates(n,1);
%[init_original main_original connector_original outit_original gate_names]=Split_Trot_Pauli_Ladder_Gates(n);
[init main connector outit gate_names]=Split_Strang_Pauli_Ladder_Gates(n);

Gates2Circ2Preview(gate,elem_gates,other_gates,['test_ladder_circ_exp_gate_n_' num2str(n)],'gate=expansion',1,150);
% Gates2Circ2Preview(gate2,elem_gates,other_gates,['test_sym_ladder_circ_exp_gate_n_' num2str(n)],'gate=expansion',1,150);
% Gates2Circ2Preview(init_original,elem_gates,other_gates,['test_ladder_circ_exp_init_orig_n_' num2str(n)],'gate=expansion',1,150);
% Gates2Circ2Preview(main_original,elem_gates,other_gates,['test_ladder_circ_exp_main_orig_n_' num2str(n)],'gate=expansion',1,150);
% Gates2Circ2Preview(connector_original,elem_gates,other_gates,['test_ladder_circ_exp_con_orig_n_' num2str(n)],'gate=expansion',1,150);
% Gates2Circ2Preview(outit_original,elem_gates,other_gates,['test_ladder_circ_exp_outit_orig_n_' num2str(n)],'gate=expansion',1,150);
Gates2Circ2Preview(init,elem_gates,other_gates,['test_ladder_circ_exp_init_n_' num2str(n)],'gate=expansion',1,150);
Gates2Circ2Preview(main,elem_gates,other_gates,['test_ladder_circ_exp_main_n_' num2str(n)],'gate=expansion',1,150);
Gates2Circ2Preview(connector,elem_gates,other_gates,['test_ladder_circ_exp_connector_n_' num2str(n)],'gate=expansion',1,150);
Gates2Circ2Preview(outit,elem_gates,other_gates,['test_ladder_circ_exp_outit_n_' num2str(n)],'gate=expansion',1,150);