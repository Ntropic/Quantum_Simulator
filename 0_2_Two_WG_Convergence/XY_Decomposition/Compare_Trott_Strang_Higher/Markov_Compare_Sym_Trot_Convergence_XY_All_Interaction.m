%Markov_Compare_Photon_Numbers_Circ_XY_All_Interaction.m
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
n2=2:6

%Trotter Iteration number
steps=8;

dt=0.0001;
E=0.005; 

%Iteration time
t=pi/4;
h=figure('Position', [10 10 1000 700]);
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
name=['comp_sym_trot_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_dt_' num2str(dt) '_E_' num2str(E) '_XY'];
sym_name=['markov_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t)  '_dt_' num2str(dt) '_E_' num2str(E) '_XY'];
trot_name=['markov_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_dt_' num2str(dt) '_E_' num2str(E)  '_XY'];
higher_name=['higher_compare_n_' num2str(min(n_h)) '_to_' num2str(max(n_h)) '_s_' num2str(s_h) '_r_' num2str(r_h) '_t_' num2str(t) '_XY'];

fprintf('Check if files already exist\n')

if exist(['../Symmetric_Trotter/' sym_name '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name '.mat']);
    F_av2_ind=a.F_av_ind;
    gate_num2=a.gate_num;
    s_q2=a.s_q;
else
    error('Files (sym Trotter) not existent')
end

%% -------------------------------------------------------------------------------------------------------------------------

if exist(['../Simple_Trotter/' trot_name '.mat'])
    a=load(['../Simple_Trotter/' trot_name '.mat']);
    F_av_ind=a.F_av_ind;
    gate_num=a.gate_num;
    s_q=a.s_q;
    fprintf('Found Simple Trotter calculations\n')
else
    error('Files (simple Trotter) non existant')
end


%% Analyze -------------------------------------------------------------------------------------------------------------


%Plot 1-Fidelity over Steps and Photons
cmap=autumn(255);
cmap2=winter(255);
P=repmat((2:max(n_m)-1)',1,steps);
S=repmat(1:steps,max(n_m)-2,1);
f_1=F_av_ind(2:end,:)';
f_2=F_av2_ind(2:end,:)';
s1=surf(ax,P',S',f_1);
hold on;
s2=surf(ax,P',S',f_2);
hold off;
xlabel('# photons')
ylabel('# Trotter steps')
zlabel('Error Rate')
legend({'Markov Simple Trotter','Markov Symmetric Trotter'})
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
set(ax,'ZScale', 'lin')
drawnow;


%Plot 1-Fidelity over gate number and Photons
h2=figure('Position', [10 10 1000 700]);
ax2=axes(h2);

cmap=autumn(255);
cmap2=winter(255);
P=repmat((2:max(n_m)-1)',1,steps);
gate_number=[gate_num(2:end,:),gate_num2(2:end,:)];
f_1=F_av_ind(2:end,:)';
f_2=F_av2_ind(2:end,:)';
s1=surf(ax2,P',gate_num(2:end,:)',f_1);
hold on;
s2=surf(ax2,P',gate_num2(2:end,:)',f_2);
hold off;
xlabel('# photons')
ylabel('# gates')
zlabel('Error Rate')
legend({'Markov Simple Trotter','Markov Symmetric Trotter'})
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
set(ax2,'ZScale', 'lin')
drawnow;
% matlab2tikz('markov_xy_trotter_comp.tex','figurehandle',h)
% matlab2tikz('markov_xy_gate_comp.tex','figurehandle',h2)

% A=getframe(h);
% B=getframe(h2);
% imwrite(A.cdata,'markov_xy_trotter_comp.png');
% imwrite(B.cdata,'markov_xy_gate_comp.png');