%Compare_Photon_Numbers_Symmetric_Pauli_Gray_All_Interaction.m
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
n_m=1:5;

%Trotter Iteration number
steps=5;

%Iteration time
t=pi/4;

name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
F_av_ind=zeros(length(n_m),steps);
gate_num=zeros(length(n_m),steps);
s_q=zeros(length(n_m),steps,2);
time_l=zeros(length(n_m),steps);
for i=1:length(n_m)
    fprintf([num2str(i) ' Photons\n'])
    n=n_m(i);
    timer=tic();
    %% Exact decomposition
    [Hij,s]=Gray_Exchange_Hamiltonian_Particles(n); %s is qubits per mode 2s is total qubit number
    %FockPrint(Hij)
    U_exact=expm(-1i*Hij*t);

    %% Approximate decomposition
    [init main connector outit gate_names]=Split_Strang_Pauli_Gray_Gates(n);
    
    other_gates2=Comp_Gate_Merger(other_gates,init,main,connector,outit);
    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    sizle=init.size;
    anc_size=init.anc_size;
    for l=1:steps
        fprintf(['  ' num2str(l) ' Trotter steps\n'])
        timer=tic();

        trott_gate=Trotter_InOut(t,l,'','',sizle,anc_size,gate_names);
        %Gates2Circ2Preview(trott_gate,elem_gates,other_gates,'test_trott_gate','gate=expansion',2,30);

        other_gates3=Comp_Gate_Merger(other_gates2,trott_gate);
        other_gates3=Gates_Tables_Prepare(other_gates3,elem_gates);   

        ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates3);
        U_approx=Expansion2Matrix(other_gates3(ind(2)));
        %average_fidelity
        indexes=Gray_Indexes(n);
        [F_av(i,l) F_av_ind(i,l)]=Average_Fidelity(U_approx,U_exact,indexes);
        time_l(n,l)=toc(timer);
        
        gate_num(i,l)=length(other_gates3(ind(2)).exp_steps.index);
        s_q(i,l,1:2)=[0,0];
        for j=1:gate_num(i,l)
            lens=length(other_gates3(ind(2)).exp_steps.index{j});
            s_q(i,l,lens)=s_q(i,lens)+1;
        end
        time(i,l)=toc(timer);
        surf(F_av_ind');
        xlabel('Number of photons')
        ylabel('Trotter steps')
        zlabel('Average Fidelity')
        view(45,45)
        drawnow;
    end 
end

matlab2tikz([name '_fidelity_surf.tex']);

save([name '.mat'],'gate_num','s_q','F_av_ind','time');

surf(gate_num');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('Number of Gates')
view(-45,45)
drawnow;

matlab2tikz([name '_gate_number_surf.tex']);