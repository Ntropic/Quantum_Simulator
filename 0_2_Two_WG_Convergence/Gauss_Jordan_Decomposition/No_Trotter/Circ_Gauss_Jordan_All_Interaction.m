%Circ_Gauss_Jordan_All_Interaction.m
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

%Photon number
n=3;

%Iteration time
t=pi/4;

%% Exact decomposition
[Hij,s]=Gray_Exchange_Hamiltonian_Particles(n); %s is qubits per mode 2s is total qubit number
%FockPrint(Hij)
a=rand(8,1);
U_exact=expm(-1i*Hij*t);
%FockPrint(U_exact)

%% Approximate decomposition

sizle=s;
if s==2
    anc_size=0;
else
    anc_size=s-2;
end
trott_gate=Fast_Gauss_Jordan_Decomposition( U_exact,1);
other_gates2=Comp_Gate_Merger(comp_gates,trott_gate);
%Gates2Circ2Preview(trott_gate,elem_gates,other_gates2,'test_trott_gate','gate=expansion',1,30);        
%U_approx=Gate2Matrix(elem_gates,comp_gates,trott_gate)
FockPrint(U_exact)

other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
U_approx=full(Expansion2Matrix(other_gates2(ind(2))));
FockPrint(U_approx)

%average_fidelity
indexes=Gray_Indexes(n);
[F_av F_av_ind]=Average_Fidelity(U_approx,U_exact,indexes)