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
n=4;

N=3;

%Trotter Iteration number
steps=5;

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
    
    U_2=expm(-1i*H*t/l);
    gj_dec_gate=Fast_Gauss_Jordan_Decomposition( U_2);
    other_gates=Comp_Gate_Merger(comp_gates,gj_dec_gate);
    %Gates2Circ2Preview(trott_gate,elem_gates,other_gates2,'test_trott_gate','gate=expansion',1,30);        

    sizle=gj_dec_gate.size;
    anc_size=gj_dec_gate.anc_size;

    trott_gate=MWG_Trotter(gj_dec_gate.names,N,l,'','',sizle,anc_size,t);
    
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