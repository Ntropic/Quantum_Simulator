%Circ_Pauli_Gray_All_Interaction.m
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

%Photon number
n=2;

N=3;

%Trotter Iteration number
steps=3;

rec_order=0;

%Iteration time
t=pi/4;

%% Exact decomposition
[Hij]=Ladder_Exchange_Hamiltonian(n);
Hij=H2MWG(Hij,N);
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);

%% Approximate decomposition
[init outit init2 outit2 symit2 con gate_names2]=MWG_Split_Strang_Pauli_Ladder_Gates(n);

other_gates=Comp_Gate_Merger(other_gates,init,outit,init2,outit2,symit2,con);
other_gates=Gates_Tables_Prepare(other_gates,elem_gates);   

sizle=2*(n+1);
anc_size=0;

connections=[(1:N-1)',(2:N)'];
weights=ones(N-1,1);
for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();
    
    pattern=rec_order;
    gate_name=['H_{L,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
    circ_str=['H_{L,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];

    trott_gate=MWG_Higher_Trotter_InOut(t,l,pattern,connections,weights,gate_name,circ_str,sizle,anc_size,gate_names2);

    other_gates2=Comp_Gate_Merger(other_gates,trott_gate);
    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   
    
    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
    U_approx=Expansion2Matrix(other_gates2(ind(2)));

    %average_fidelity
    indexes=Ladder_Indexes(n,N);
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,indexes);
    time_l(l)=toc(timer);
end 

f1_cmp_ind=1-F_av_ind;

figure()
%subplot(3,1,1:2);
plot(1:steps,f1_cmp_ind,'r')
%xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^2)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;

%subplot(3,1,3);
%plot(1:steps,time_l,'r')
%xlabel('Trotter steps')
%ylabel('t [s]')

%drawnow;

figure()
n=2:steps;
logn=log((n-1)./n);
u1_diff=log(f1_cmp_ind(2:end))-log(f1_cmp_ind(1:end-1));
semilogx(n,u1_diff./logn,'r')
xlabel('Trotter steps')
ylabel('Convergence Order')
legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence_order1_50.tex')
drawnow;