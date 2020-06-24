%Compare_Photon_Numbers_Circ_Gauss_Jordan_All_Interaction.m
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
n_m=1:7;

%Iteration time
t=pi/4;

for i=1:length(n_m)
    n=n_m(i);
    %% Exact decomposition
    timer=tic();
    [Hij,s]=Gray_Exchange_Hamiltonian_Particles(n); %s is qubits per mode 2s is total qubit number
    %FockPrint(Hij)
    a=rand(8,1);
    %a([1,3,4,5,6,8])=zeros(6,1);
    %U_exact=diag(exp(1i*a));
    U_exact=expm(-1i*Hij*t);
    %FockPrint(U_exact)

    %% Approximate decomposition

    sizle=s;
    if s==2
        anc_size=0;
    else
        anc_size=s-2;
    end
    trott_gate=Fast_Gauss_Jordan_Decomposition( U_exact);
    other_gates2=Comp_Gate_Merger(comp_gates,trott_gate);
    %Gates2Circ2Preview(trott_gate,elem_gates,other_gates2,'test_trott_gate','gate=expansion',1,30);        
    %U_approx=Gate2Matrix(elem_gates,comp_gates,trott_gate)
    %FockPrint(U_exact)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
    U_approx=full(Expansion2Matrix(other_gates2(ind(2)),[],1));
    %FockPrint(U_approx)
    gate_num(i)=length(other_gates2(ind(2)).exp_steps.index);
    s_q(i,1:2)=[0,0];
    for j=1:gate_num(i)
        lens=length(other_gates2(ind(2)).exp_steps.index{j});
        s_q(i,lens)=s_q(i,lens)+1;
    end
    
    %average_fidelity
    indexes=Gray_Indexes(n);
    [F_av(i) F_av_ind(i)]=Average_Fidelity(U_approx,U_exact,indexes);
    time(i)=toc(timer);
    
    
    %subplot(3,1,1);
    %plot(n_m(1:i),F_av(1:i));
    %xlabel('F_av');
    %subplot(3,1,2);
    %plot(n_m(1:i),time(1:i));
    %xlabel('t [s]');
    %subplot(3,1,3);
    area(n_m(1:i),s_q(1:i,1)+s_q(1:i,2),'FaceColor','red');
    hold on;
    area(n_m(1:i),s_q(1:i,1),'FaceColor','blue');
    legend({'# two qubit gates','# single qubit gates'})
    xlabel('# of photons in modes');
    hold off;
    drawnow;
    
    %figure()
    %plot(n_m(1:i),s_q(1:i,1)./gate_num(1:i)')
end

save('compare_n_1_to_7_gauss_jordan.mat','gate_num','s_q','F_av_ind','time');

matlab2tikz('gate_numbers_abs.tex');

set(gca, 'YScale', 'log');

matlab2tikz('gate_numbers_abs_log.tex');   
    
%--------------------------------------------------------------------------

figure()
plot(n_m(1:i),s_q(1:i,1)./gate_num(1:i)')
xlabel('# of photons in modes');
ylabel('# 1 qubit gates / # 2 qubit gates')

matlab2tikz('gate_numbers_ratio.tex');