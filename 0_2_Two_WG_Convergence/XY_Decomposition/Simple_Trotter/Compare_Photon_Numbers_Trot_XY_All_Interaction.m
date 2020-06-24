%Compare_Photon_Numbers_Trot_XY_All_Interaction.m
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
n_m=2:10;

%Trotter Iteration number
steps=20;

%Iteration time
t=pi/4;

h=figure();
ax=axes(h);

name_files=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];

if exist([name_files '.mat'])==0
    F_av_ind=zeros(length(n_m),steps);
    gate_num=zeros(length(n_m),steps);
    s_q=zeros(length(n_m),steps,2);
    time_l=zeros(length(n_m),steps);
    for i=1:length(n_m)
        n=n_m(i);
        N=n-1;
        fprintf([num2str(N) ' Photons\n'])

        timer=tic();
        %% Exact decomposition
        [Hij]=XY_Exchange_Hamiltonian(n);
        %FockPrint(Hij)
        U_exact=expm(-1i*Hij*t);

        %% Prepare gates
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
        %A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

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
        %B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);


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
            [gate,step_gate trotter_name trotter_step_name]=Trotter(t,l,name,circ_string,circ_string_step,n,0,gates_names,gates_indexes);
            %Do Trotter steps
            other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
            other_gates=Gates_Tables_Prepare(other_gates,elem_gates);
            ind=Gate_Index_by_Name(trotter_name,elem_gates,other_gates);
            U_approx=Expansion2Matrix(other_gates(ind(2)));

            %average_fidelity
            indexes=2.^(0:n-1)+1;
            [F_av(i,l) F_av_ind(i,l)]=Average_Fidelity(U_approx,U_exact,indexes);
            time_l(n,l)=toc(timer);

            gate_num(i,l)=length(other_gates(ind(2)).exp_steps.index);
            s_q(i,l,1:2)=[0,0];
            for j=1:gate_num(i,l)
                lens=length(other_gates(ind(2)).exp_steps.index{j});
                s_q(i,l,lens)=s_q(i,l,lens)+1;
            end
            time(i,l)=toc(timer);
            surf(ax,F_av_ind');
            xlabel('Number of photons')
            ylabel('Trotter steps')
            zlabel('Average Fidelity')
            view(45,45)
            drawnow;
        end 
    end
    save([name_files '.mat'],'gate_num','s_q','F_av_ind','time');
else
    fprintf('Already calculated. Loading data\n')
    load([name_files '.mat'],'gate_num','s_q','F_av_ind','time');
end
       

P=repmat((2:max(n_m)-1)',1,steps);
S=repmat(1:steps,max(n_m)-2,1);
surf(ax,P',S',1-F_av_ind(2:end,:)');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('1-Average Fidelity')
set(gca,'ZScale', 'log')
set(gca,'xlim',[min(P(:)),max(P(:))]);
set(gca,'ylim',[min(S(:)),max(S(:))]);
view(45,-45)
drawnow;
        
matlab2tikz([name_files '_1-fidelity_surf.tex']);

h2=figure();
ax2=axes(h2);
Pfull=repmat((1:max(n_m)-1)',1,steps);
Sfull=repmat(1:steps,max(n_m)-1,1);
surf(ax2,Pfull',Sfull',gate_num');
xlabel('Number of photons')
ylabel('Trotter steps')
zlabel('Number of Gates')
set(gca,'xlim',[min(Pfull(:)),max(Pfull(:))]);
set(gca,'ylim',[min(Sfull(:)),max(Sfull(:))]);
view(-45,45)
drawnow;

matlab2tikz([name_files '_gate_number_surf.tex']);