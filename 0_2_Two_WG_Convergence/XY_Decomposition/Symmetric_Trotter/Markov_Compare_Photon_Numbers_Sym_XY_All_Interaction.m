%Compare_Photon_Numbers_Sym_XY_All_Interaction.m
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

%Qubit number (0ne bigger than photon number)
n_m=2:9;

%Trotter Iteration number
steps=8;

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

h=figure();
ax=axes(h);

name_files=['markov_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_dt_' num2str(t_step) '_E_' num2str(E) '_XY'];
if exist([name_files '.mat'])==0
    F_av_ind=zeros(length(n_m),steps);
    gate_num=zeros(length(n_m),steps);
    s_q=zeros(length(n_m),steps,2);
    time_l=zeros(length(n_m),steps);
    for i=1:length(n_m)
        n=n_m(i);
        N=n-1;
        fprintf([num2str(N) ' Photons\n'])
        
        %Create initial pure state
        m=n;
        theta=0;    %relative phase
        phi=NOON_Wave(1,n,m,theta);
        rho0=full(sparse(Wave2Density(phi,1)));
        
        timer=tic();
        %% Exact decomposition
        [Hij]=XY_Exchange_Hamiltonian(n);
        %FockPrint(Hij)
        U_exact=expm(-1i*Hij*t);

        %Subspace for analysis of algorithm
        indexes=2.^(0:n-1)+1;

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

            %% Approximate decomposition
            %Create Trotter steps
            name=['XY_n_' num2str(n)];
            circ_string=['XY_{\text{approx}}^{n=' num2str(n) '}'];
            circ_string_step=['XY_{\text{step}}^{n=' num2str(n) '}'];
            gates_names={nameA,nameB};
            gates_indexes={[1:A.size],1+[1:B.size]};
            other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);

            %This command creates a Trotter gate for all the steps -> Great for
            %creating a gate circuit
            [gate,step_gate,trotter_name,trotter_step_name]=Strang_Trotter(t,l,name,circ_string,circ_string_step,n,0,gates_names,gates_indexes);
            %if l==2
            %    Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',1,30);
            %end
            %Do Trotter steps
            other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
            other_gates=Gates_Tables_Prepare(other_gates,elem_gates);
            ind=Gate_Index_by_Name(trotter_name,elem_gates,other_gates);
            [rho sizes anc_sizes]=Noisy_Expansion2Matrix(rho0,errors,other_gates(ind(2)));
            rho_i=U_exact*rho0*U_exact';

            %average_fidelity
            [F_av(i,l) F_av_ind(i,l)]=Hilbert_Schmidt(rho_i,rho,indexes);
            time_l(n,l)=toc(timer);

            gate_num(i,l)=length(other_gates(ind(2)).exp_steps.index);
            s_q(i,l,1:2)=[0,0];
            for j=1:gate_num(i,l)
                lens=length(other_gates(ind(2)).exp_steps.index{j});
                s_q(i,l,lens)=s_q(i,lens)+1;
            end
            time(i,l)=toc(timer);
            
            surf(ax,F_av_ind');
            xlabel('Number of photons')
            ylabel('Trotter steps')
            zlabel('Average error')
            view(45,45)
            drawnow;
        end 
    end
    save([name_files '.mat'],'gate_num','s_q','F_av_ind','time');
else
    fprintf('Already calculated. Loading data\n')
    load([name_files '.mat'],'gate_num','s_q','F_av_ind','time');
end

P=repmat((2:max(n_m)-1)',1,steps)-1;
S=repmat(1:steps,max(n_m)-2,1);
surf(ax,P',S',1-F_av_ind(2:end,:)');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('1-Error rate')
set(gca,'ZScale', 'log')
set(gca,'xlim',[min(P(:)),max(P(:))]);
set(gca,'ylim',[min(S(:)),max(S(:))]);
view(45,22.5)
drawnow;
        
matlab2tikz([name_files '_1-error_rate_surf.tex']);