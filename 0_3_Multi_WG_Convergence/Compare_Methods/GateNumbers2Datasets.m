%Waveguides_Increase_All_GateNumbers2Datasets.m
%How many more Trotter steps are needed when adding another Waveguide
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

for j=1:2
    fprintf([num2str(j) 'th dataset: \n'])
    if j==1
        %Photon number
        n_m=1:7;

        %Trotter Iteration number
        steps=10;

        N_m=2:4;
    elseif j==2
        %Photon number
        n_m=1:3;

        %Trotter Iteration number
        steps=10;

        N_m=4:6;
    end

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

        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat']);
            F_trot=a.F_av_ind;
            fprintf(' - Found Simple Trotter calculations\n')
        else
            error('Not found (trott fid)')
        end
        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat']);
            s_q_trot=a.s_q;
            fprintf(' - Found Simple Trotter (gate number) calculations\n')
        else
            error('Not found (trott gate numbers)')
        end

        %% -------------------------------------------------------------------------------------------------------------------------
        if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits{i} '.mat']);
            F_sym=a.F_av_ind;
            fprintf(' - Found Simple Trotter calculations\n')
        else
            error('Not found (sym fid)')
        end
        if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat']);
            s_q_sym=a.s_q;
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

            if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits{i} '.mat'])
                a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits{i} '.mat']);
                F_high_r0{counter}=a.F_av_ind;
                fprintf([' - Found Higher Order Trotter [r=' num2str(r) '] calculations\n'])
            else
                error('Not found (higher fid)')
            end
            if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat'])
                a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat']);
                s_q_high_r0{counter}=a.s_q;
                fprintf([' - Found Higher_Order_Trotter [r=' num2str(r) '] (gate number) calculations\n'])
            else
                error(['Not found (Higher_Order_Trotter) [r=' num2str(r) ']'])
            end
        end
        
        if j==1
            if i==1
                F={};
                s_q={};
            end
            F={F{:},F_trot,F_sym,F_high_r0{:}};
            s_q={s_q{:},s_q_trot,s_q_sym,s_q_high_r0{:}};
            if i==1
                len_sq=length(s_q);
            end
        else
            if i==1
                F2={};
                s_q2={};
            end
            F2={F2{:},F_trot,F_sym,F_high_r0{:}};
            s_q2={s_q2{:},s_q_trot,s_q_sym,s_q_high_r0{:}};
        end
    end
end

F=reshape(F,len_sq,length(exits));
F=F';
s_q=reshape(s_q,len_sq,length(exits));
s_q=s_q';

F2=reshape(F2,len_sq,length(exits));
F2=F2';
s_q2=reshape(s_q2,len_sq,length(exits));
s_q2=s_q2';

fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%% Analyze Data---------------------------------------------------------------------------------------------------------

for method=1:3
    h_c_wg=figure('Position',[100,100,1400,600]);

    %Get data
    s_q_ind=s_q{method,1};
    s_q2_ind=s_q2{method,1};

    s_q_ind=s_q_ind(:,:,1,1)+s_q_ind(:,:,1,2);
    s_q2_ind=s_q2_ind(:,:,1,1)+s_q2_ind(:,:,1,2);


    surf([1:7],[2:4],s_q_ind);
    hold on;
    surf([1:3],[4:6],s_q2_ind(3:end,:));

    xlabel('$N_\text{max}$')
    zlabel({'$n(N_\text{max})$'})
    ylabel('$M$')

    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[1:7];
    set(gca, 'ZScale', 'log');
    view(-45,60)
    rotate3d on;
    axis tight;
    drawnow;
    matlab2tikz(['gate_numbers_MN_' exits{method} '.tex'],'standalone',true,'parsestrings',false,'width','6cm','height','4cm')
end