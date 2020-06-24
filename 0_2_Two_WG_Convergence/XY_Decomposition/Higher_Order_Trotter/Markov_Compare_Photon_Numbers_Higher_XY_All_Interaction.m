%Markov_Compare_Photon_Numbers_Circ_Pauli_Gray_All_Interaction.m
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
n_m=1:6;

%Trotter Iteration number
steps=5;

%Higher Trotter recursion order
rec_order=0;

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

name=['markov_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_r_' num2str(rec_order) '_t_' num2str(t) '_dt_' num2str(t_step) '_E_' num2str(E) '_XY'];
F_av_ind=zeros(length(n_m)-1,steps);
gate_num=zeros(length(n_m)-1,steps);
s_q=zeros(length(n_m)-1,steps,2);
time_l=zeros(length(n_m)-1,steps);
for i=2:length(n_m)
    n=n_m(i);
    N=n-1;%Photon number
    fprintf([num2str(N) ' Photons\n'])
    timer=tic();
          
    %Create initial pure state
    m=n;
    theta=0;    %relative phase
    phi=NOON_Wave(1,n,m,theta);
    rho0=full(sparse(Wave2Density(phi,1)));
        
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

        pattern=ones(1,rec_order);
        gate_name=['H_{XY,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
        circ_string=['H_{XY,n=' num2str(n) '}^{(' num2str(rec_order+1) ')}'];
        gates_names={nameA,nameB};
        gates_indexes={[1:A.size],1+[1:B.size]};
        other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
        %This command creates a Trotter gate for all the steps -> Great for
        %creating a gate circuit
        if rec_order==1
            pattern=ones(1,1);
        else
            pattern=zeros(1,1);
        end
        [seq,d_seq,gate,step_gate]=Higher_Trotter(t,l,pattern,'',gate_name,circ_string,n,anc_size,gates_names,gates_indexes);
        %if l==1
        %    Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',1,30);
        %end

        %Do Trotter steps
        other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
        other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

        ind=Gate_Index_by_Name(gate.names,elem_gates,other_gates);
        [rho sizes anc_sizes]=Noisy_Expansion2Matrix(rho0,errors,other_gates(ind(2)));
        rho_i=U_exact*rho0*U_exact';
        %average_fidelity
        indexes=2.^(0:n-1)+1;
        [F_av(i-1,l) F_av_ind(i-1,l)]=Hilbert_Schmidt(rho_i,rho,indexes);
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
        zlabel('Average Error')
        view(45,45)
        drawnow;
    end 
end

matlab2tikz([name '_error_surf.tex']);

save([name '.mat'],'gate_num','s_q','F_av_ind','time');

P=repmat((2:max(n_m))',1,steps)-1;
S=repmat(1:steps,max(n_m)-1,1);
surf(P',S',1-F_av_ind(1:end,:)');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('1-Error rate')
set(gca,'ZScale', 'log')
set(gca,'xlim',[min(P(:)),max(P(:))]);
set(gca,'ylim',[min(S(:)),max(S(:))]);
view(45,-45)
drawnow;
        
% matlab2tikz([name_files '_1-error_rate_surf.tex']);