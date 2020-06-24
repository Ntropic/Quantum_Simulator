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
load('../XY_quantum_gates.mat','-mat')

%Photon number
n_m=2:6;

%Trotter Iteration number
steps=15;

%Higher Trotter recursion order
rec_order=[3];
rec_len=1;

%Iteration time
t=pi/4;

if rec_len==1
    name=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_r_' num2str(rec_order,'%1i_') 't_' num2str(t) '_XY'];
else
    name=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_r_' num2str(rec_order,'%1i_') 'r_l_' num2str(rec_len,'%1i_') 't_' num2str(t) '_XY'];
end
F_av_ind=zeros(length(n_m)-1,steps);
gate_num=zeros(length(n_m)-1,steps);
s_q=zeros(length(n_m)-1,steps,2);
time_l=zeros(length(n_m)-1,steps);
for i=length(n_m)-1:length(n_m)
    n=n_m(i);
    N=n-1;%Photon number
    fprintf([num2str(N) ' Photons\n'])
    timer=tic();
    %% Exact decomposition
    
    Hij=XY_Exchange_Hamiltonian(n);
    U_exact=expm(-1i*t*Hij);

    %Subspace for analysis of algorithm
    index=2.^(0:n-1)+1;

    Ni=N:-1:1;
    in=1:N;
    Jx=sqrt(in).*sqrt(Ni);

    %Gate A
    sizle=floor(n/2)*2;
    anc_size=0;
    nameA=['A_n_' num2str(n)];
    circ_string=['A_{XY}^{n=' num2str(n) '}'];
    gates_indexes={};
    param_steps={};
    gates_names={};
    for k=1:2:2*floor(n/2)-1
        gates_indexes={gates_indexes{:},[k,k+1]};
        param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k)}]};
    end
    [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
    A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
    A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

    %Gate B
    sizle=floor((n-1)/2)*2;
    anc_size=0;
    nameB=['B_n_' num2str(n)];
    circ_string=['B_{XY}^{n=' num2str(n) '}'];
    gates_indexes={};
    param_steps={};
    gates_names={};
    for k=1:2:2*floor((n-1)/2)-1
        gates_indexes={gates_indexes{:},[k,k+1]};
        param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k+1)}]};
    end
    [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
    B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
    B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);   

    for l=1:steps
        fprintf(['  ' num2str(l) ' Trotter steps\n'])
        timer=tic();
        
        if length(rec_order)>1
            pattern=rec_order;
        else
            pattern=rec_order*ones(1,rec_len);
        end
        gate_name=['H_{XY,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
        circ_string=['H_{XY,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
        gates_names={nameA,nameB};
        gates_indexes={[1:A.size],1+[1:B.size]};
        other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
        %This command creates a Trotter gate for all the steps -> Great for
        %creating a gate circuit
        [seq,d_seq,gate,step_gate]=Higher_Trotter(t,l,pattern,'',gate_name,circ_string,n,anc_size,gates_names,gates_indexes);
        %if l==1
        %    Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',1,30);
        %end

        %Do Trotter steps
        other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
        other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

        ind=Gate_Index_by_Name(gate.names,elem_gates,other_gates);
        U_approx=Expansion2Matrix(other_gates(ind(2)));
        %average_fidelity
        indexes=2.^(0:n-1)+1;
        [F_av(i-1,l) F_av_ind(i-1,l)]=Average_Fidelity(U_approx,U_exact,indexes);
        time_l(n,l)=toc(timer);
        
        gate_num(i-1,l)=length(other_gates(ind(2)).exp_steps.index);
        s_q(i,l,1:2)=[0,0];
        for j=1:gate_num(i-1,l)
            lens=length(other_gates(ind(2)).exp_steps.index{j});
            s_q(i,l,lens)=s_q(i,lens)+1;
        end
        time(i-1,l)=toc(timer);
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