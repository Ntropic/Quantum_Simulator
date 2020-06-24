clear all;
close all;
clc;
%Plot the errors for different methods after 1 trotter step with 
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

    r=0:3;

    %Iteration time
    t=pi/4;

    %% Fidelities and gate numbers (gn)
    %% ------------------------------------------------------------------------------------------------------------------------
    exits='_pauli_gray';
    exit_Folders={'Pauli_Ladder_Decomposition','Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};

    trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    fprintf('-> Check if simulation results are located in the designated folders:\n\n')

    %% -------------------------------------------------------------------------------------------------------------------------
        if exist(['../Simple_Trotter/' trot_name exits '.mat'])
            a=load(['../Simple_Trotter/' trot_name exits '.mat']);
            F_trot=a.F_av_ind;
            fprintf(' - Found Simple Trotter calculations\n')
        else
            error('Not found (trott)')
        end
        if exist(['../Simple_Trotter/' trot_name_gn exits '.mat'])
            a=load(['../Simple_Trotter/' trot_name_gn exits '.mat']);
            s_q_trot=a.s_q;
            fprintf(' - Found Simple Trotter (gate number) calculations\n')
        else
            error('Not found (trott gate numbers)')
        end

        %% -------------------------------------------------------------------------------------------------------------------------
        if exist(['../Symmetric_Trotter/' sym_name exits '.mat'])
            a=load(['../Symmetric_Trotter/' sym_name exits '.mat']);
            F_sym=a.F_av_ind;
            fprintf(' - Found Symetric Trotter calculations\n')
        else
            error('Not found (sym)');
        end
        if exist(['../Symmetric_Trotter/' sym_name_gn exits '.mat'])
            a=load(['../Symmetric_Trotter/' sym_name_gn exits '.mat']);
            s_q_sym=a.s_q;
            fprintf(' - Found Symetric Trotter (gate number) calculations\n')
        else
            error('Not found (sym gate numbers)');
        end

        %% -------------------------------------------------------------------------------------------------------------------------
        c=0;
        for r_h=r
           	higher_name_r0=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h) '_N_' num2str(max(N_m))];
            higher_name_r0_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h) '_N_' num2str(max(N_m))];

            c=c+1;
            if exist(['../Higher_Order_Trotter/' higher_name_r0 exits '.mat'])
                a=load(['../Higher_Order_Trotter/' higher_name_r0 exits '.mat']);
                F_high_r0{c}=a.F_av_ind;
                fprintf(' - Found Higher_Order_Trotter [r0] calculations\n')
            else
                error('Not found (Higher_Order_Trotter) [r0]')
            end
            if exist(['../Higher_Order_Trotter/' higher_name_r0_gn exits '.mat'])
                a=load(['../Higher_Order_Trotter/' higher_name_r0_gn exits '.mat']);
                s_q_high_r0{c}=a.s_q;
                fprintf(' - Found Higher_Order_Trotter [r0] (gate number) calculations\n')
            else
                error('Not found (Higher_Order_Trotter) [r0]')
            end
        end
        %% -------------------------------------------------------------------------------------------------------------------------

    if j==1
        F={F_trot,F_sym,F_high_r0{:}};
        s_q={s_q_trot,s_q_sym,s_q_high_r0{:}};
    elseif j==2
        %First shorten second dataset
        F2={};
        s_q2={};
        F2={F_trot,F_sym,F_high_r0{:}};
        s_q2={s_q_trot,s_q_sym,s_q_high_r0{:}};
        for i=1:length(F2)
            A=F2{i};
            A=A(3:5,:,:);
            F2{i}=A;
            B=s_q2{i};
            B=B(3:5,:,:,:);
            s_q2{i}=B;    
        end
        n_m2=1:3;
        N_m2=4:6;
        n_m=1:7;
        N_m=2:4;
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

%% Errors for different particle numbers
step=5;
names={'$\hat{T}$','$\hat{S}_1$','$\hat{S}_3^{(3)}$','$\hat{S}_3^{(5)}$','$\hat{S}_3^{(7)}$','$\hat{S}_3^{(9)}$'};
%stri={'-.r','r','--k','--b','--m','--c'};
%stri2={':r',':r',':k',':b',':m',':c'};
len_sq=length(s_q);
for i=1:len_sq
    %First dataset
    f=F{i};
    s=s_q{i};
    t=s(: ,: ,step,1)+s(:,:,step,2);
    g=f(:,:,step);
    g=1-g;
    g(1,1)=0;
    
    %Second Dataset
    f2=F2{i};
    s2=s_q2{i};
    t2=s2(: ,: ,step,1)+s2(:,:,step,2);
    g2=f2(:,:,step);
    g2=1-g2;
    
    
    figure()
    surf(n_m,N_m,g)
    hold on;
    surf(n_m2,N_m2,g2)
    set(gca, 'ZScale', 'log');
    rotate3d on;
    xlabel('\# Particles N')
    ylabel('\# of Modes M')
    zlabel('$1-F_\text{av}$')
    xlim([1 7])
    ylim([2 6])
    hold off;
end
%legend(names,'interpreter','latex')
%matlab2tikz(['pauli_gray_errors_N_' num2str(N) '.tex'],'standalone',true,'parseStrings',false)
