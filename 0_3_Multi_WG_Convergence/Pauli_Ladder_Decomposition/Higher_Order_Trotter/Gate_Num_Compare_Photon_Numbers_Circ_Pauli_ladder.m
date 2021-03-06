%Gate_Num_Compare_Photon_Numbers_Circ_Pauli_ladder.m
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

%Trotter Iteration number
steps=10;

for op=3;

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
    %Higher Trotter recursion order
    for rec_order=0:3;

        %Iteration time
        t=pi/4;

        name=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(rec_order) '_N_' num2str(max(N_m)) '_pauli_ladder'];
        F_av_ind=zeros(length(n_m),steps);
        gate_num=zeros(length(n_m),steps);
        s_q=zeros(length(N_m),length(n_m),steps,2);
        time_l=zeros(length(n_m),steps);
        for h=1:length(N_m)
            N=N_m(h);
            fprintf([num2str(N) ' Waveguides\n'])

            connections=[(1:N-1)',(2:N)'];
            weights=ones(N-1,1);
            for i=1:length(n_m)
                fprintf(['  ' num2str(i) ' Photons\n'])
                n=n_m(i);
                %% Exact decomposition
                %[Hij]=Ladder_Exchange_Hamiltonian(n);
                %FockPrint(Hij)
                %U_exact=expm(-1i*Hij*t);

                %% Approximate decomposition
                [init outit init2 outit2 symit con gate_names]=MWG_Split_Strang_Pauli_Ladder_Gates(n);

                other_gates2=Comp_Gate_Merger(other_gates,init,outit,init2,outit2,symit,con);

                sizle=init.size;
                anc_size=init.anc_size;
                for l=1:steps
                    fprintf(['    ' num2str(l) ' Trotter steps\n'])

                    pattern=rec_order;
                    gate_name=['H_{L,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
                    circ_str=['H_{L,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];

                    trott_gate=MWG_Higher_Trotter_InOut(t,l,pattern,connections,weights,gate_name,circ_str,sizle,anc_size,gate_names);
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

        save([name '.mat'],'gate_num','s_q');
    end
end