%mwg_gates_creator_test.m
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
load('cxy_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

%-------------------------------------------------------------------------
n=1;
M=9;
steps=2;

%gate=Pauli_Gray_Gates(n);
[init,outit,init2,outit2,main,connector,gate_names]=MWG_Split_Strang_Pauli_Gray_Gates(n);
other_gates=Comp_Gate_Merger(other_gates,init,outit,init2,outit2,main,connector);

gate=MWG_Trotter_InOut(pi/4,M,1,steps,[],[],init.size,init.anc_size,gate_names)
%gate=MWG_Strang_Trotter_InOut(pi/4,M,1,steps,[],[],init.size,init.anc_size,gate_names)
Gates2Circ2Preview(gate,elem_gates,other_gates,['mwg_higher_test_M_' num2str(M) '_' num2str(n)],'gate=expansion',1,150);

%Gates2Circ2Preview(gate,elem_gates,other_gates,['test_gray_circ_exp_gate_n_' num2str(n)],'gate=expansion',1,150);

%Gates2Circ2Preview(init,elem_gates,other_gates,['test_gray_circ_exp_init_n_' num2str(n)],'gate=expansion',1,150);
%Gates2Circ2Preview(main,elem_gates,other_gates,['test_gray_circ_exp_main_n_' num2str(n)],'gate=expansion',1,150);
%Gates2Circ2Preview(connector,elem_gates,other_gates,['test_gray_circ_exp_connector_n_' num2str(n)],'gate=expansion',1,150);
%Gates2Circ2Preview(outit,elem_gates,other_gates,['test_gray_circ_exp_outit_n_' num2str(n)],'gate=expansion',1,150);
%Gates2Circ2Preview(init2,elem_gates,other_gates,['test_gray_circ_exp_init2_n_' num2str(n)],'gate=expansion',1,150);
%Gates2Circ2Preview(outit2,elem_gates,other_gates,['test_gray_circ_exp_outit2_n_' num2str(n)],'gate=expansion',1,150);