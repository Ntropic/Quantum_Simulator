%gates_sep_test.m
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
load('../../../../Gates_Table/elem_gates.mat','-mat')
load('../../../../Gates_Table/comp_gates.mat','-mat')

%Photon number
n=1;

%Iteration time
t=pi/4;

%% Exact decomposition
[Hij,s]=Gray_Exchange_Hamiltonian_Particles(n); %s is qubits per mode 2s is total qubit number
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);

%% Approximate decomposition

sizle=s;
if s==2
    anc_size=0;
else
    anc_size=s-2;
end
trott_gate=Fast_Gauss_Jordan_Decomposition( U_exact);
trott_gate.steps
fprintf('-------------------------\n')
for i=1:trott_gate.step_num
    trott_gate.steps.gates{i}
    [trott_gate.steps.param{i}{:}]
    trott_gate2=Create_Comp_Gate('testgate',2,0,{trott_gate.steps.gates{i}},{trott_gate.steps.index{i}},{trott_gate.steps.param{i}})
    other_gates2=Comp_Gate_Merger(comp_gates,trott_gate2);
    %Gates2Circ2Preview(trott_gate,elem_gates,other_gates2,'test_trott_gate','gate=expansion',2,30);        
    U_approx=Gate2Matrix(elem_gates,comp_gates,trott_gate2)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(trott_gate2.names,elem_gates,other_gates2);
    U_approx2=full(Expansion2Matrix(other_gates2(ind(2))))
    
    fprintf('--------------------------------\n')
end