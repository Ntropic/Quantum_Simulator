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

%Photon number
n=7;

N=3;

%Trotter Iteration number
steps=5;

rec_order=0;
%Iteration time
t=pi/4;

%% Exact decomposition
[H]=Gray_Exchange_Hamiltonian_Particles(n);
Hij=H2MWG(H,N);   %Create Multi Waveguide Gate
    
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);

%% Exact decomposition

for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();
    
    [other_gates2 sizle anc_size names names2 names_sum]=MWG_Gauss_Jordan_Decomposition_Prepare(H,t,l,rec_order,comp_gates);

    trott_gate=MWG_Higher_Trotter(names2,names,names_sum,N,l,'','',sizle,anc_size,t);
    
    other_gates2=Comp_Gate_Merger(other_gates2,trott_gate);
    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   
    
    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
    U_approx=Expansion2Matrix(other_gates2(ind(2)));
    

    %average_fidelity
    indexes=Gray_Indexes(n,N);
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,indexes);
    time_l(l)=toc(timer);
end 

f1_cmp_ind=1-F_av_ind;

figure();
%subplot(3,1,1:2);
plot(1:steps,f1_cmp_ind,'r')
%xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^2)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;