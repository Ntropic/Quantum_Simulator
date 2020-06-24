%Symmetric_Circ_Pauli_Gray_All_Interaction.m
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
load('clifford_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

%Photon number
n=2;

%Trotter Iteration number
steps=5;

N=3;

%Iteration time
t=pi/2;

%% Exact decomposition
[Hij]=Ladder_Exchange_Hamiltonian(n);
Hij=H2MWG(Hij,3);
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);

%% Approximate decomposition
[init outit init2 outit2 symit con gate_names2]=MWG_Split_Strang_Pauli_Ladder_Gates(n);


other_gates_mwg=Comp_Gate_Merger(other_gates,init,outit,init2,outit2,symit,con);
other_gates_mwg=Gates_Tables_Prepare(other_gates_mwg,elem_gates);   


sizle=2*(n+1);
anc_size=0;
for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();
    
    trott_gate=MWG_Strang_Trotter_InOut(t,N,[],l,'','',sizle,anc_size,gate_names2 );
    %Gates2Circ2Preview(trott_gate,elem_gates,other_gates,'test_trott_gate','gate=expansion',2,30);

    
    other_gates_2=Comp_Gate_Merger(other_gates_mwg,trott_gate);
    other_gates_2=Gates_Tables_Prepare(other_gates_2,elem_gates);   
    
    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates_2);
    U_approx=Expansion2Matrix(other_gates_2(ind(2)));

    %average_fidelity
    indexes=Ladder_Indexes(n,N);
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,indexes);
    
    time_l(l)=toc(timer);
end 

f1_cmp=1-F_av_ind

figure()
subplot(3,1,1:2);
plot(1:steps,f1_cmp,'r')
hold on;
plot(1:steps,f2_cmp,'r')
%xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^4)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;
