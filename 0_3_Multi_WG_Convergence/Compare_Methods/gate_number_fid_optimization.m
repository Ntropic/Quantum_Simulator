%gate_number_fid_optimization.m
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

%Convergence Plot
N=6;
M=4;


%Photon number
n_m=1:7;
N_m=2:4;

modes_max=4;
particles_max=7;


%Trotter Iteration number
steps=10;

r_h=0:3;



%Iteration time
t=pi/4;

%% Fidelities and gate numbers (gn)
%% ------------------------------------------------------------------------------------------------------------------------
exits2={'_pauli_gray','_gj'};
exit_Folders={'Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};
for i=1:2
    exits=exits2{i};
    
    trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    fprintf('-> Check if simulation results are located in the designated folders:\n\n')

    %% -------------------------------------------------------------------------------------------------------------------------
    if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits '.mat'])
        a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits '.mat']);
        F_trot=a.F_av_ind;
        fprintf(' - Found Simple Trotter calculations\n')
    else
        error('Not found (trott calc)')
    end
    if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits '.mat'])
        a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits '.mat']);
        s_q_trot=a.s_q;
        fprintf(' - Found Simple Trotter (gate number) calculations\n')
    else
        error('Not found (trott gate numbers)')
    end

    %% -------------------------------------------------------------------------------------------------------------------------
    if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits '.mat'])
        a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits '.mat']);
        F_sym=a.F_av_ind;
        fprintf(' - Found Symetric Trotter calculations\n')
    else
        error('Not found (sym gcalc)');
    end
    if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits '.mat'])
        a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits '.mat']);
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

        if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits '.mat'])
            a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits '.mat']);
            F_high_r0{counter}=a.F_av_ind;
            fprintf(' - Found Higher_Order_Trotter calculations\n')
        else
            error(['Not found (Higher_Order_Trotter calc) [r=' num2str(r) ']'])
        end
        if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits '.mat'])
            a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits '.mat']);
            s_q_high_r0{counter}=a.s_q;
            fprintf([' - Found Higher_Order_Trotter [r=' num2str(r) '] (gate number) calculations\n'])
        else
            error(['Not found (Higher_Order_Trotter gate numbers) [r=' num2str(r) ']'])
        end
    end
    if i==1
       	F={F_trot,F_sym,F_high_r0{:}};
        s_q={s_q_trot,s_q_sym,s_q_high_r0{:}};
        len_sq=length(s_q);
    else
        F2={F_trot,F_sym,F_high_r0{:}};
        s_q2={s_q_trot,s_q_sym,s_q_high_r0{:}};
    end
end


    
fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%each element of f=F{i} contains a 3 dimensional tensor with the
%dimensions f(a,b,c) are 
%       a=number of modes N_m
%       b=number of particles n_m
%       c=number of trotter steps

%each element of s=s_q{i} contains a 4 dimensional tensor with the
%dimensions s(a,b,c,d) are 
%       a=number of modes N_m
%       b=number of particles n_m
%       c=number of trotter steps
%       d==1 -> single qubit gates
%       d==2 -> two qubit gates


%% Errors over gates

for particles=2:particles_max
    for modes=3:modes_max
        figure()
        %% Gray coding
        all_F=[];
        all_s=[];
        for i=1:length(s_q)
            curr_F=F{i};
            c_F=curr_F(modes-1,particles,:);
            curr_s=s_q{i};
            c_s=curr_s(modes-1,particles,:,1)+curr_s(modes-1,particles,:,2);

            all_F=[all_F;1-c_F(:)];
            all_s=[all_s;c_s(:)];
        end
        [all_s, order]=sort(all_s);
        all_F=all_F(order);

        i=2;
        min_F=all_F(1);
        while i<length(all_F)
            if all_F(i)>min_F
                all_F(i)=[];
                all_s(i)=[];
            else
                min_F=all_F(i);
                i=i+1;
            end
        end
        semilogy(all_s,all_F,'b')
        hold on;

        %% GJ
        all_F=[];
        all_s=[];
        for i=1:length(s_q)
            curr_F=F2{i};
            c_F=curr_F(modes-1,particles,:);
            curr_s=s_q2{i};
            c_s=curr_s(modes-1,particles,:,1)+curr_s(modes-1,particles,:,2);

            all_F=[all_F;1-c_F(:)];
            all_s=[all_s;c_s(:)];
        end
        [all_s, order]=sort(all_s);
        all_F=all_F(order);

        i=2;
        min_F=all_F(1);
        while i<length(all_F)
            if all_F(i)>min_F
                all_F(i)=[];
                all_s(i)=[];
            else
                min_F=all_F(i);
                i=i+1;
            end
        end
        semilogy(all_s,all_F,'r')
        title(['$M=' num2str(modes) ', N_\text{max}=' num2str(particles) '$']);
    end
end


% xlim(x_lim);
% ylim(y_lim);
% ylabel('$1-F$');
% xlabel('\# Gates');
% legend(names,'interpreter','latex')
% matlab2tikz(['gj_errors_N_' num2str(N) '_M_' num2str(M) '.tex'],'standalone',true,'parseStrings',false)