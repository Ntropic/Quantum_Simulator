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


%Trotter Iteration number
steps=10;

for op=3
    if op==1
        n_m=1:7;
        N_m=2:4;
    elseif op==2
        n_m=1:3;
        N_m=2:6;
    else
        n_m=1:7;
        N_m=2:6;
    end

    for rec_order=0:3;

        %Iteration time
        t=pi/4;

        name=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(rec_order) '_N_' num2str(max(N_m)) '_gj'];
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

                timer=tic();

                [other_gates2 sizle anc_size names names2 names_sum]=MWG_Gauss_Jordan_Decomposition_Prepare(H,t,l,rec_order,comp_gates);

                for h=1:length(N_m)
                    N=N_m(h);
                    fprintf(['    ' num2str(N) ' Waveguides\n'])

                    connections=[(1:N-1)',(2:N)'];
                    weights=ones(N-1,1);

                    trott_gate=MWG_Higher_Trotter(names2,names,names_sum,N,l,'','',sizle,anc_size,t);

                    other_gates2=Comp_Gate_Merger(other_gates2,trott_gate);

                    ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
                    [elem_gates,other_gates2]=Gate_Counter(elem_gates, other_gates2);

                    gate=other_gates2(ind(2));
                    s_q(h,i,l,1)=gate.single_gates;
                    s_q(h,i,l,2)=gate.double_gates;

                end
            end 
        end

        save([name '.mat'],'s_q');
    end
end