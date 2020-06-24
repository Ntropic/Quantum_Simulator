%Compare_Photon_Numbers_Circ_Pauli_Gray_All_Interaction.m
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
n_m=1:4;

%Trotter Iteration number
steps=5;

%Iteration time
t=pi/4;
%Error model
t_step=0.0001;
T1=1;
T2=1;

phi=0;
E=0.005;    %Error rate
delta=sqrt(10/3*E); 
E1=5/4*E;
errors=[t_step,T1,T2,phi,delta,E1];

name=['markov_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_dt_' num2str(t_step) '_E_' num2str(E) '_pauli_ladder'];
F_av_ind=zeros(length(n_m),steps);
gate_num=zeros(length(n_m),steps);
s_q=zeros(length(n_m),steps,2);
time_l=zeros(length(n_m),steps);
for i=1:length(n_m)
    fprintf([num2str(i) ' Photons\n'])
    n=n_m(i);
    timer=tic();
                
    %Create initial pure state (in pauli_ladder_coding)
    m=n;
    theta=0;    %relative phase
    phi=NOON_Wave_Ladder(n,2,m,theta);
    rho0=full(sparse(Wave2Density(phi,1)));
    
    %% Exact decomposition
    [Hij]=Ladder_Exchange_Hamiltonian(n);
    %FockPrint(Hij)
    U_exact=expm(-1i*Hij*t);

    %% Approximate decomposition
    [init main connector outit gate_names]=Split_Strang_Pauli_Ladder_Gates(n);

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
        [rho sizes anc_sizes]=Noisy_Expansion2Matrix(rho0,errors,other_gates3(ind(2)));
        rho_i=U_exact*rho0*U_exact';
        %average_fidelity
        indexes=Ladder_Indexes(n);
        [F_av(i,l) F_av_ind(i,l)]=Hilbert_Schmidt(rho_i,rho,indexes);
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
        zlabel('Average Error')
        view(45,45)
        drawnow;
    end 
end

matlab2tikz([name '_fidelity_surf.tex']);

save([name '.mat'],'gate_num','s_q','F_av_ind','time');


P=repmat((2:max(n_m)-1)',1,steps)-1;
S=repmat(1:steps,max(n_m)-2,1);
surf(ax,P',S',1-F_av_ind(2:end,:)');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('1-Error rate')
set(gca,'ZScale', 'log')
set(gca,'xlim',[min(P(:)),max(P(:))]);
set(gca,'ylim',[min(S(:)),max(S(:))]);
view(45,-45)
drawnow;
        
matlab2tikz([name_files '_1-error_rate_surf.tex']);