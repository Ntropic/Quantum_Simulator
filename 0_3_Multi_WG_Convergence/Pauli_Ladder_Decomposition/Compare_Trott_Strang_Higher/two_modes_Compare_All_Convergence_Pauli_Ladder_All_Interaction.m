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

%Photon number
n_m=1:7;

%Trotter Iteration number
steps=10;

N_m=2:6;

r_h=0:3;

%Iteration time
t=pi/4;

%% Fidelities and gate numbers (gn)
%% ------------------------------------------------------------------------------------------------------------------------
exits={'_pauli_ladder','_pauli_gray','_gj'};
exit_Folders={'Pauli_Ladder_Decomposition','Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};

trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

fprintf('-> Check if simulation results are located in the designated folders:\n\n')

%% -------------------------------------------------------------------------------------------------------------------------
for i=1:length(exits) %Only search for gate numbers
    fprintf(['- Searching for calculations of type ' exit_Folders{i} '\n'])
    if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat'])
        a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat']);
        s_q_trot{i}=a.s_q;
        fprintf(' - Found Simple Trotter (gate number) calculations\n')
    else
        error('Not found (trott gate numbers)')
    end

    %% -------------------------------------------------------------------------------------------------------------------------
    if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat'])
        a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat']);
        s_q_sym{i}=a.s_q;
        fprintf(' - Found Symetric Trotter (gate number) calculations\n')
    else
        error('Not found (sym gate numbers)');
    end

    %% -------------------------------------------------------------------------------------------------------------------------
    counter=0;
    for r=r_h
        counter=counter+1;
        higher_name_r0=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];
        higher_name_r0_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];

        if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat']);
            s_q_high_r0{counter,i}=a.s_q;
            fprintf([' - Found Higher_Order_Trotter [r=' num2str(r) '] (gate number) calculations\n'])
        else
            error(['Not found (Higher_Order_Trotter) [r=' num2str(r) ']'])
        end
    end

end

s_q={};
for i=1:length(exits)
    s_q={s_q{:},s_q_trot{i},s_q_sym{i},s_q_high_r0{:,i}};
    if i==1
        len_sq=length(s_q);
    end
end

s_q=reshape(s_q,len_sq,length(exits));
s_q=s_q';
    
fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')
