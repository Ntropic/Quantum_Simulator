%quantum_gates.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))


%Load gates
load('../../Gates_Table/elem_gates.mat','-mat');
load('../../Gates_Table/comp_gates.mat','-mat');

seq=[1,0,0,1];

%newer_gate.size=3;
%newer_gate.anc_size=0;
%newer_gate.step_num=1;
%newer_gate.steps.index={[1,2,3]};
%newer_gate.steps.gates={'toffoli00'};
%newer_gate.steps.param={[]};

%[matrix,matrix_anc]=Gate2Matrix(elem_gates,comp_gates,newer_gate,2);
%[M index M_mode]=Sparse2Square(matrix_anc,{'fock unit',1,3});
%fprintf(M_mode)

new_gate=Cn_U2(elem_gates,comp_gates,seq,[0,pi/4,pi/4,pi/4]);
fprintf('Created gate structure\n')
new_gate.steps.gates{:}
new_gate.steps.index{:}
new_gate.steps.param{:}

[matrix,matrix_anc]=Gate2Matrix(elem_gates,comp_gates,new_gate,0);

%fprintf('With ancilla qbits:\n')
%[M_anc index_anc M_mode_anc]=Sparse2Square(matrix_anc,{'fock unit',1,4,2});
%fprintf(M_mode_anc)

fprintf('\nWithout ancilla qbits:\n')
[M index M_mode]=Sparse2Square(matrix,'fock unit');
fprintf(M_mode)