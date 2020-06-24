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

N_m=2:4;

r_h=0;
r_h2=1;

%Iteration time
t=pi/4;

%% Fidelities and gate numbers (gn)
%% ------------------------------------------------------------------------------------------------------------------------
trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_ladder'];
trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_ladder'];

sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_ladder'];
sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_ladder'];

higher_name_r0=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h) '_N_' num2str(max(N_m)) '_pauli_ladder'];
higher_name_r0_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h) '_N_' num2str(max(N_m)) '_pauli_ladder'];

higher_name_r1=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h2) '_N_' num2str(max(N_m)) '_pauli_ladder'];
higher_name_r1_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r_h2) '_N_' num2str(max(N_m)) '_pauli_ladder'];


fprintf('-> Check if simulation results are located in the designated folders:\n\n')

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Simple_Trotter/' trot_name '.mat'])
    a=load(['../Simple_Trotter/' trot_name '.mat']);
    F_trot=a.F_av_ind;
    fprintf(' - Found Simple Trotter calculations\n')
else
    error('Not found (trott)')
end
if exist(['../Simple_Trotter/' trot_name_gn '.mat'])
    a=load(['../Simple_Trotter/' trot_name_gn '.mat']);
    s_q_trot=a.s_q;
    fprintf(' - Found Simple Trotter (gate number) calculations\n')
else
    error('Not found (trott gate numbers)')
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Symmetric_Trotter/' sym_name '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name '.mat']);
    F_sym=a.F_av_ind;
    fprintf(' - Found Symetric Trotter calculations\n')
else
    error('Not found (sym)');
end
if exist(['../Symmetric_Trotter/' sym_name_gn '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name_gn '.mat']);
    s_q_sym=a.s_q;
    fprintf(' - Found Symetric Trotter (gate number) calculations\n')
else
    error('Not found (sym gate numbers)');
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Higher_Order_Trotter/' higher_name_r0 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r0 '.mat']);
    F_high_r0=a.F_av_ind;
    fprintf(' - Found Higher_Order_Trotter [r0] calculations\n')
else
    error('Not found (Higher_Order_Trotter) [r0]')
end
if exist(['../Higher_Order_Trotter/' higher_name_r0_gn '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r0_gn '.mat']);
    s_q_high_r0=a.s_q;
    fprintf(' - Found Higher_Order_Trotter [r0] (gate number) calculations\n')
else
    error('Not found (Higher_Order_Trotter) [r0]')
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Higher_Order_Trotter/' higher_name_r1 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r1 '.mat']);
    F_high_r1=a.F_av_ind;
    fprintf(' - Found Higher_Order_Trotter [r1] calculations\n')
else
    error('Not found (Higher_Order_Trotter) [r1]')
end
if exist(['../Higher_Order_Trotter/' higher_name_r1_gn '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r1_gn '.mat']);
    s_q_high_r1=a.s_q;
    fprintf(' - Found Higher_Order_Trotter [r1] (gate number) calculations\n')
else
    error('Not found (Higher_Order_Trotter) [r1]')
end

fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%% Analyze Data---------------------------------------------------------------------------------------------------------

%% Fastest Fidelity
F={F_trot,F_sym,F_high_r0,F_high_r1};
s_q={s_q_trot,s_q_sym,s_q_high_r0,s_q_high_r1};
names={'Trotter','Symmetric','Tripling','Quintupling'};

err=0.1;   %Desired error rate (upper bound)
fid=1-err; %Desired fidelity (lower bound)
costs=[0.01,1]; %[single gate cost,double gate cost]

%Prepare Predictors
fprintf(' - Predicting Gate numbers beyond simulation data. \n')
s_q_pred=Difference_Predictor(s_q,3);
fprintf(' - Predicting Gate Fidelities beyond simulation data. \n')
F_pred=Order_Predictor(F);
%Get Fastest Fidelity
fprintf(' - Createing map of costs to achieve a desired fidelity. \n')
cost_maps=Fastest_Fidelity(fid,F,s_q,costs,F_pred,s_q_pred);
fprintf(' - Plotting Fastest Fidelity map. \n')
mode='redraw_pixel';
data={s_q_pred,F_pred,F,s_q,costs,names,N_m,n_m};
h=Fastest_Fidelity_Plot([],mode,cost_maps,names,N_m,n_m,fid,data);



% %% Plot 1-Fidelity over Steps and  and WG's
% %Colormaps
% cmap_trot=autumn(255);
% cmap2=winter(255);
% cmap3=cool(255);
% cmap4=parula(255);
% %Error rates (1-Fidelity)
% f_trot=1-F_trot(2:end,:)';
% f_sym=1-F_sym(2:end,:)';
% f_high_r0=1-F_high_r0(2:end,:)';
% f_high_r1=1-F_high_r1(2:end,:)';
% 
% P=repmat((1:max(n_m)-1)',1,steps);
% S=repmat(1:steps,max(n_m)-1,1);
% 
% s1=surf(ax,P',S',f_trot);
% hold on;
% s2=surf(ax,P',S',f_sym);
% s3=surf(ax,P',S',f_3);
% s4=surf(ax,P',S',f_4);
% hold off;
% xlabel('# photons')
% ylabel('# Trotter steps')
% zlabel('1-Average Fidelity')
% legend({'Simple','Symmetric','Tripling','Quintupling'})
% set(ax,'ZScale', 'log')
% set(ax,'xlim',[min(P(:)),max(P(:))]);
% set(ax,'ylim',[min(S(:)),max(S(:))]);
% view(-135,45)
% limz=log10(get(ax,'zlim'));
% minz=min(limz);
% limz=limz-minz;
% f_trot=log10(f_trot)-minz;
% f_sym=log10(f_sym)-minz;
% dlimz=limz(2)/length(cmap_trot);
% f_1_i=floor(f_trot/dlimz)+1;
% f_2_i=floor(f_sym/dlimz)+1;
% f_3_i=floor(f_3/dlimz)+1;
% f_4_i=floor(f_4/dlimz)+1;
% f_1_c=reshape(cmap_trot(f_1_i(:,:),:),[size(f_trot,1) size(f_trot,2) 3]);
% f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_sym,1) size(f_sym,2) 3]);
% f_3_c=reshape(cmap3(f_3_i(:,:),:),[size(f_3,1) size(f_3,2) 3]);
% f_4_c=reshape(cmap4(f_4_i(:,:),:),[size(f_4,1) size(f_4,2) 3]);
% 
% s1.CData=f_1_c;
% s2.CData=f_2_c;
% s3.CData=f_3_c;
% s4.CData=f_4_c;
% set(gcf,'color',0.98*[1 1 1])
% A=getframe(h);
% imwrite(A.cdata,'higher2_ladder_trotter_comp.png');
% drawnow;
% 
% %Plot 1-Fidelity over gate number and Photons
% h2=figure('Position', [10 10 1000 700]);
% ax2=axes(h2);
% 
% P=repmat((1:max(n_m)-1)',1,steps);
% gate_number=[gate_num(2:end,:),gate_num2(2:end,:)];
% gate_number=[gate_number(:)',gate_num3(:)',gate_num4(:)'];
% f_trot=1-F_av_ind(2:end,:)';
% f_sym=1-F_av2_ind(2:end,:)';
% f_3=1-F_av3_ind(2:end,:)';
% f_4=1-F_av4_ind(2:end,:)';
% s1=surf(ax2,P',gate_num(2:end,:)',f_trot);
% hold on;
% s2=surf(ax2,P',gate_num2(2:end,:)',f_sym);
% s3=surf(ax2,P',gate_num3(2:end,:)',f_3);
% s4=surf(ax2,P',gate_num4(2:end,:)',f_4);
% hold off;
% xlabel('# photons')
% ylabel('# gates')
% zlabel('1-Average Fidelity')
% legend({'Simple','Symmetric','Tripling','Quintupling'})
% set(ax2,'ZScale', 'log')
% set(ax2,'YScale','log')
% set(ax2,'xlim',[min(P(:)),max(P(:))]);
% set(ax2,'ylim',[min(gate_number(:)),max(gate_number(:))]);
% view(-135,45);
% limz=log10(get(ax2,'zlim'));
% minz=min(limz);
% limz=limz-minz;
% eps=10^-6;
% f_trot=log10(eps+f_trot)-minz;
% f_sym=log10(eps+f_sym)-minz;
% f_4=log10(eps+f_4)-minz;
% f_3=log10(eps+f_3)-minz;
% dlimz=limz(2)/length(cmap_trot);
% f_1_i=floor(f_trot/dlimz)+1;
% f_2_i=floor(f_sym/dlimz)+1;
% f_3_i=floor(f_3/dlimz)+1;
% f_4_i=floor(f_4/dlimz)+1;
% f_1_c=reshape(cmap_trot(f_1_i(:,:),:),[size(f_trot,1) size(f_trot,2) 3]);
% f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_sym,1) size(f_sym,2) 3]);
% f_3_c=reshape(cmap3(f_3_i(:,:),:),[size(f_3,1) size(f_3,2) 3]);
% f_4_c=reshape(cmap4(f_4_i(:,:),:),[size(f_4,1) size(f_4,2) 3]);
% s1.CData=f_1_c;
% s2.CData=f_2_c;
% s3.CData=f_3_c;
% s4.CData=f_4_c;
% set(gcf,'color',0.98*[1 1 1])
% drawnow;
% % matlab2tikz('xy_trotter_comp.tex','figurehandle',h)
% % matlab2tikz('xy_gate_comp.tex','figurehandle',h2)
% 
% B=getframe(h2);
% imwrite(B.cdata,'higher2_ladder_gate_comp.png');