%Gates_M2_GJ_2Datasets.m
%Compare Gray and Ladder gates
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

    %Photon number
    n_m=1:18;%1:7;

    %Trotter Iteration number
    steps=1;%0;

    N_m=2;%2:4;

    r_h=0:3;

    %Iteration time
    t=pi/4;

    %% Fidelities and gate numbers (gn)
    %% ------------------------------------------------------------------------------------------------------------------------
    exits={'_pauli_ladder','_pauli_gray','_gj'};%,'_gj2'};
    exit_Folders={'Pauli_Ladder_Decomposition','Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};%,'Gauss_Jordan_Decomposition'};

    trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    fprintf('-> Check if simulation results are located in the designated folders:\n\n')

    %% -------------------------------------------------------------------------------------------------------------------------
    for i=1:length(exits) %Only search for gate numbers
        fprintf(['- Searching for calculations of type ' exit_Folders{i} '\n'])

%         if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat'])
%             a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat']);
%             F_trot=a.F_av_ind;
%             fprintf(' - Found Simple Trotter calculations\n')
%         else
%             error('Not found (trott fid)')
%         end
        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat']);
            s_q_trot=a.s_q;
            fprintf(' - Found Simple Trotter (gate number) calculations\n')
        else
            error('Not found (trott gate numbers)')
        end
        
        if i==1
%            F={};
            s_q={};
        end
%        F={F{:},F_trot};
        s_q={s_q{:},s_q_trot};
        if i==1
            len_sq=length(s_q);
        end
    end

%F=reshape(F,len_sq,length(exits));
%F=F';
s_q=reshape(s_q,len_sq,length(exits));
s_q=s_q';


fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%% Analyze Data---------------------------------------------------------------------------------------------------------
x=1:18;
%Ladder
s_q_ind=s_q{1,1}; %Ladder
s_q_1=s_q_ind(1,:,1,1);
s_q_2=s_q_ind(1,:,1,2);
s_q_2l=s_q_2+s_q_1;

h=figure('Position',[100,100,800,600]);    
fill([x x(end:-1:1)],[s_q_2+s_q_1 zeros(size(x))],[1 1 0]);
hold on;
fill([x x(end:-1:1)],[s_q_1 zeros(size(x))],[0 0 1]);

%set(gca, 'YScale', 'log')
xlabel('$N_\text{max}$')
ylabel('$n(N_\text{max})$')


s_q_ind=s_q{2,1}; %Gray
s_q_1g=s_q_ind(1,:,1,1);
s_q_2g=s_q_ind(1,:,1,2);
s_q_2g=s_q_2g+s_q_1g;

ax=plot(s_q_2g);
ax.Color=[1 0.5 0];
y_lim=ylim();
axis([1 18 y_lim]);

legend({'Two qubits','Single qubit','Gray'},'interpreter','latex','location','northwest');
matlab2tikz(['gate_number_ladder.tex'],'standalone',true,'parsestrings',false,'width','6cm','height','4cm')

%% Now Gray Plot
s_q_ind=s_q{2,1}; %Gray
s_q_1=s_q_ind(1,:,1,1);
s_q_2=s_q_ind(1,:,1,2);
s_q_2g=s_q_2+s_q_1;

h2=figure('Position',[100,100,800,600]);    
fill([x x(end:-1:1)],[s_q_2+s_q_1 zeros(size(x))],[1 1 0]);
hold on;
fill([x x(end:-1:1)],[s_q_1 zeros(size(x))],[0 0 1]);
%set(gca, 'YScale', 'log')
xlabel('$N_\text{max}$')
ylabel('$n(N_\text{max})$')
y_lim=ylim();
axis([1 18 y_lim]);

legend({'Two qubits','Single qubit'},'interpreter','latex','location','northwest');
matlab2tikz(['gate_number_gray.tex'],'standalone',true,'parsestrings',false,'width','6cm','height','4cm')


%% How many Trotter steps are necessary to get to the same amount of gates
a=s_q_2l; %Ladder
b=s_q_2g; %Gray
c=a./b;
d=c(2:18);
h3=figure('Position',[100,100,800,600]);    
plot(2:18,d,'k');
y_lim=ylim(gca);
hold on;
plot([3 3],y_lim,':k')
plot([7 7],y_lim,':k')
plot([15 15],y_lim,':k')
axis([2 18 y_lim])
xlabel('$N_\text{max}$')
ylabel('$\frac{n_{\text{L}}(N_\text{max})}{n_{\text{G}}(N_\text{max})}$')
z=2:18;

matlab2tikz(['ratio_gate_number_g_l.tex'],'standalone',true,'parsestrings',false,'width','6cm','height','4cm')