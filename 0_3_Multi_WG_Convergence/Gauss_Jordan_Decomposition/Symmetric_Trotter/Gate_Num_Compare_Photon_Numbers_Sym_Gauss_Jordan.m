%Gate_Num_Compare_Photon_Numbers_Circ_Gauss_Jordan.m
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

%Trotter Iteration number
steps=10;

N_m=2:6;

%Iteration time
t=pi/4;

name=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_gj'];
F_av_ind=zeros(length(n_m),steps);
gate_num=zeros(length(n_m),steps);
s_q=zeros(length(N_m),length(n_m),steps,2);
time_l=zeros(length(n_m),steps);

for i=1:length(n_m)
    fprintf([num2str(i) ' Photons\n'])
    n=n_m(i);

    [H]=Gray_Exchange_Hamiltonian_Particles(n);
    %% Exact decomposition
    %Hij=H2MWG(H,N);
    %FockPrint(Hij)
    %U_exact=expm(-1i*Hij*t);

    for l=1:steps
        fprintf(['  ' num2str(l) ' Trotter steps\n'])

        name2_gate=['U_n_' num2str(n) '_t_' num2str(t/l)];
        circ_name2=['U_{n=' num2str(n) '}(t=' num2str(t/l) ')'];
        U_2=expm(-1i*H*t/l);
        gj_dec_gate2=Fast_Gauss_Jordan_Decomposition(U_2,0,name2_gate);

        name_gate=['U_n_' num2str(n) '_t_' num2str(t/l/2)];
        circ_name2=['U_{n=' num2str(n) '}(t=' num2str(t/l/2) ')'];
        U=expm(-1i*H*t/l/2);              
        gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U,0,name_gate);

        other_gates2=Comp_Gate_Merger(comp_gates,gj_dec_gate,gj_dec_gate2);

        sizle=gj_dec_gate.size;
        anc_size=gj_dec_gate.anc_size;

        for h=1:length(N_m)
            N=N_m(h);
            fprintf(['    ' num2str(N) ' Waveguides\n'])

            connections=[(1:N-1)',(2:N)'];
            weights=ones(N-1,1);

            trott_gate=MWG_Strang_Trotter({name_gate,name2_gate},N,l,'','',sizle,anc_size,t);
            %Gates2Circ2Preview(trott_gate,elem_gates,other_gates,'test_trott_gate','gate=expansion',2,30);

            other_gates3=Comp_Gate_Merger(other_gates2,trott_gate); 

            ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates3);
            [elem_gates,other_gates3]=Gate_Counter(elem_gates, other_gates3);

            gate=other_gates3(ind(2));
            s_q(h,i,l,1)=gate.single_gates;
            s_q(h,i,l,2)=gate.double_gates;

        end
    end 
end

save([name '.mat'],'s_q');
