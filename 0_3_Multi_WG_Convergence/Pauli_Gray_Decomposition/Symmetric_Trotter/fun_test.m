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
load('../cxy_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

%Photon number
n=2;

N=3;

%Trotter Iteration number
steps=5;

%Iteration time
t=pi/4;

%% Exact decomposition
[Hij]=Gray_Exchange_Hamiltonian_Particles(n);
Hij=H2MWG(Hij,3);   %Create Multi Waveguide Gate
    
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);


%% Approximate decomposition
[init outit init2 outit2 symit con gate_names]=MWG_Split_Strang_Pauli_Gray_Gates(n);

other_gates=Comp_Gate_Merger(other_gates,init,outit,init2,outit2,symit,con);
other_gates=Gates_Tables_Prepare(other_gates,elem_gates);   

sizle=init.size;
anc_size=init.anc_size;
for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();
    
    trott_gate=MWG_Strang_Trotter_InOut(t,N,[],l,'','',sizle,anc_size,gate_names);
    
    other_gates2=Comp_Gate_Merger(other_gates,trott_gate);
    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   
    
    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
    U_approx=Expansion2Matrix(other_gates2(ind(2)));
    

    %average_fidelity
    indexes=Gray_Indexes(n,N);
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