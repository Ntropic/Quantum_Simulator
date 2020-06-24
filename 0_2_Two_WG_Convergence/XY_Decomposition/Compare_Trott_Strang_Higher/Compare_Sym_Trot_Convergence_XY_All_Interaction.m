%Compare_Photon_Numbers_Circ_XY_All_Interaction.m
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
steps=10;

%Iteration time
t=pi/4;
h=figure('Position', [10 10 1000 700]);
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
name=['comp_sym_trot_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];
sym_name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];
trot_name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];

fprintf('Check if files already exist\n')

if exist(['../Symmetric_Trotter/' sym_name '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name '.mat']);
    F_av2_ind=a.F_av_ind;
    gate_num2=a.gate_num;
    s_q2=a.s_q;
else
    fprintf('Calculate symmetric Trotter calculations\n')
    F_av2_ind=zeros(length(n_m),steps);
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
            U_approx=Expansion2Matrix(other_gates(ind(2)));

            %average_fidelity
            [F_av(i,l) F_av_ind(i,l)]=Average_Fidelity(U_approx,U_exact,indexes);
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
            zlabel('Average Fidelity')
            view(45,45)
            drawnow;
        end 
    end
    save(['../Symmetric_Trotter/' sym_name '.mat'],'gate_num','s_q','F_av_ind','time');
    F_av2_ind=F_av_ind;
    s_q2=s_q;
    gate_num2=gate_num;
    time2=time;
end

%% -------------------------------------------------------------------------------------------------------------------------

if exist(['../Simple_Trotter/' trot_name '.mat'])
    a=load(['../Simple_Trotter/' trot_name '.mat']);
    F_av_ind=a.F_av_ind;
    gate_num=a.gate_num;
    s_q=a.s_q;
    fprintf('Found Simple Trotter calculations\n')
else
    fprintf('Calculate simple Trotter calculations\n')
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
                s_q(i,l,lens)=s_q(i,lens)+1;
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
    save(['../Simple_Trotter/' trot_name '.mat'],'gate_num','s_q','F_av_ind','time');
end


%% Analyze -------------------------------------------------------------------------------------------------------------


%Plot 1-Fidelity over Steps and Photons
cmap=autumn(255);
cmap2=winter(255);
P=repmat((2:max(n_m)-1)',1,steps);
S=repmat(1:steps,max(n_m)-2,1);
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
s1=surf(ax,P',S',f_1);
hold on;
s2=surf(ax,P',S',f_2);
hold off;
xlabel('# photons')
ylabel('# Trotter steps')
zlabel('1-Average Fidelity')
legend({'Simple Trotter','Symmetric Trotter'})
set(ax,'ZScale', 'log')
set(ax,'xlim',[min(P(:)),max(P(:))]);
set(ax,'ylim',[min(S(:)),max(S(:))]);
view(-135,45)
limz=log10(get(ax,'zlim'));
minz=min(limz);
limz=limz-minz;
f_1=log10(f_1)-minz;
f_2=log10(f_2)-minz;
dlimz=limz(2)/length(cmap);
f_1_i=floor(f_1/dlimz)+1;
f_2_i=floor(f_2/dlimz)+1;
f_1_c=reshape(cmap(f_1_i(:,:),:),[size(f_1,1) size(f_1,2) 3]);
f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_2,1) size(f_2,2) 3]);
s1.CData=f_1_c;
s2.CData=f_2_c;
set(gcf,'color',0.98*[1 1 1])
drawnow;


%Plot 1-Fidelity over gate number and Photons
h2=figure('Position', [10 10 1000 700]);
ax2=axes(h2);

cmap=autumn(255);
cmap2=winter(255);
P=repmat((2:max(n_m)-1)',1,steps);
gate_number=[gate_num(2:end,:),gate_num2(2:end,:)];
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
s1=surf(ax2,P',gate_num(2:end,:)',f_1);
hold on;
s2=surf(ax2,P',gate_num2(2:end,:)',f_2);
hold off;
xlabel('# photons')
ylabel('# gates')
zlabel('1-Average Fidelity')
legend({'Simple Trotter','Symmetric Trotter'})
set(ax2,'ZScale', 'log')
set(ax2,'YScale','log')
set(ax2,'xlim',[min(P(:)),max(P(:))]);
set(ax2,'ylim',[min(gate_number(:)),max(gate_number(:))]);
view(-135,45);
limz=log10(get(ax2,'zlim'));
minz=min(limz);
limz=limz-minz;
f_1=log10(f_1)-minz;
f_2=log10(f_2)-minz;
dlimz=limz(2)/length(cmap);
f_1_i=floor(f_1/dlimz)+1;
f_2_i=floor(f_2/dlimz)+1;
f_1_c=reshape(cmap(f_1_i(:,:),:),[size(f_1,1) size(f_1,2) 3]);
f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_2,1) size(f_2,2) 3]);
s1.CData=f_1_c;
s2.CData=f_2_c;
set(gcf,'color',0.98*[1 1 1])
drawnow;
% matlab2tikz('xy_trotter_comp.tex','figurehandle',h)
% matlab2tikz('xy_gate_comp.tex','figurehandle',h2)

A=getframe(h);
B=getframe(h2);
imwrite(A.cdata,'xy_trotter_comp.png');
imwrite(B.cdata,'xy_gate_comp.png');