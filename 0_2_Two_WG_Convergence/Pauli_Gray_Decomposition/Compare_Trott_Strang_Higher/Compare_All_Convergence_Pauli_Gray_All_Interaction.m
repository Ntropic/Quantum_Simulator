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
n_m=1:5;
n_h=1:4;

%Trotter Iteration number
steps=5;
r_h=0;
r_h2=1;

%Iteration time
t=pi/4;

h=figure('Position', [10 10 1000 700]);
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
name=['comp_sym_trot_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
sym_name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
trot_name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
higher_name=['higher_compare_n_' num2str(min(n_h)) '_to_' num2str(max(n_h)) '_s_' num2str(steps) '_r_' num2str(r_h) '_t_' num2str(t) '_pauli_gray'];
higher_name2=['higher_compare_n_' num2str(min(n_h)) '_to_' num2str(max(n_h)) '_s_' num2str(steps) '_r_' num2str(r_h2) '_t_' num2str(t) '_pauli_gray'];
fprintf('Check if files already exist\n')

if exist(['../Symmetric_Trotter/' sym_name '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name '.mat']);
    F_av2_ind=a.F_av_ind;
    gate_num2=a.gate_num;
    s_q2=a.s_q;
else
    error('Not found (sym)');
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Simple_Trotter/' trot_name '.mat'])
    a=load(['../Simple_Trotter/' trot_name '.mat']);
    F_av_ind=a.F_av_ind;
    gate_num=a.gate_num;
    s_q=a.s_q;
    fprintf('Found Simple Trotter calculations\n')
else
    error('Not found (trott)')
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Higher_Order_Trotter/' higher_name '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name '.mat']);
    F_av3_ind=a.F_av_ind;
    gate_num3=a.gate_num;
    s_q3=a.s_q;
    fprintf('Found Higher_Order_Trotter calculations\n')
else
    error('Not found (Higher_Order_Trotter)')
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Higher_Order_Trotter/' higher_name2 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name2 '.mat']);
    F_av4_ind=a.F_av_ind;
    gate_num4=a.gate_num;
    s_q4=a.s_q;
    fprintf('Found Higher_Order_Trotter2 calculations\n')
else
    error('Not found (Higher_Order_Trotter2)')
end


%% Analyze -------------------------------------------------------------------------------------------------------------


%Plot 1-Fidelity over Steps and Photons
cmap=autumn(255);
cmap2=winter(255);
cmap3=cool(255);
cmap4=parula(255);
P=repmat((1:max(n_m)-1)',1,steps);
S=repmat(1:steps,max(n_m)-1,1);
P2=repmat((1:max(n_h)-1)',1,steps);
S2=repmat(1:steps,max(n_h)-1,1);
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
f_3=1-F_av3_ind(2:end,:)';
f_4=1-F_av4_ind(2:end,:)';
s1=surf(ax,P',S',f_1);
hold on;
s2=surf(ax,P',S',f_2);
s3=surf(ax,P2',S2',f_3);
s4=surf(ax,P2',S2',f_4);
hold off;
xlabel('# photons')
ylabel('# Trotter steps')
zlabel('1-Average Fidelity')
legend({'Simple','Symmetric','Tripling','Quintupling'})
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
f_3_i=floor(f_3/dlimz)+1;
f_4_i=floor(f_4/dlimz)+1;
f_1_c=reshape(cmap(f_1_i(:,:),:),[size(f_1,1) size(f_1,2) 3]);
f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_2,1) size(f_2,2) 3]);
f_3_c=reshape(cmap3(f_3_i(:,:),:),[size(f_3,1) size(f_3,2) 3]);
f_4_c=reshape(cmap4(f_4_i(:,:),:),[size(f_4,1) size(f_4,2) 3]);

s1.CData=f_1_c;
s2.CData=f_2_c;
s3.CData=f_3_c;
s4.CData=f_4_c;
set(gcf,'color',0.98*[1 1 1])
A=getframe(h);
imwrite(A.cdata,'higher2_pauli_gray_trotter_comp.png');
drawnow;

%Plot 1-Fidelity over gate number and Photons
h2=figure('Position', [10 10 1000 700]);
ax2=axes(h2);

P=repmat((1:max(n_m)-1)',1,steps);
P2=repmat((1:max(n_h)-1)',1,steps);
gate_number=[gate_num(2:end,:),gate_num2(2:end,:)];
gate_number=[gate_number(:)',gate_num3(:)',gate_num4(:)'];
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
f_3=1-F_av3_ind(2:end,:)';
f_4=1-F_av4_ind(2:end,:)';
s1=surf(ax2,P',gate_num(2:end,:)',f_1);
hold on;
s2=surf(ax2,P',gate_num2(2:end,:)',f_2);
s3=surf(ax2,P2',gate_num3(2:end,:)',f_3);
s4=surf(ax2,P2',gate_num4(2:end,:)',f_4);
hold off;
xlabel('# photons')
ylabel('# gates')
zlabel('1-Average Fidelity')
legend({'Simple','Symmetric','Tripling','Quintupling'})
set(ax2,'ZScale', 'log')
set(ax2,'YScale','log')
set(ax2,'xlim',[min(P(:)),max(P(:))]);
set(ax2,'ylim',[min(gate_number(:)),max(gate_number(:))]);
view(-135,45);
limz=log10(get(ax2,'zlim'));
minz=min(limz);
limz=limz-minz;
eps=10^-6;
f_1=log10(eps+f_1)-minz;
f_2=log10(eps+f_2)-minz;
f_4=log10(eps+f_4)-minz;
f_3=log10(eps+f_3)-minz;
dlimz=limz(2)/length(cmap);
f_1_i=floor(f_1/dlimz)+1;
f_2_i=floor(f_2/dlimz)+1;
f_3_i=floor(f_3/dlimz)+1;
f_4_i=floor(f_4/dlimz)+1;
f_1_c=reshape(cmap(f_1_i(:,:),:),[size(f_1,1) size(f_1,2) 3]);
f_2_c=reshape(cmap2(f_2_i(:,:),:),[size(f_2,1) size(f_2,2) 3]);
f_3_c=reshape(cmap3(f_3_i(:,:),:),[size(f_3,1) size(f_3,2) 3]);
f_4_c=reshape(cmap4(f_4_i(:,:),:),[size(f_4,1) size(f_4,2) 3]);
s1.CData=f_1_c;
s2.CData=f_2_c;
s3.CData=f_3_c;
s4.CData=f_4_c;
set(gcf,'color',0.98*[1 1 1])
drawnow;
% matlab2tikz('xy_trotter_comp.tex','figurehandle',h)
% matlab2tikz('xy_gate_comp.tex','figurehandle',h2)

B=getframe(h2);
imwrite(B.cdata,'higher2_gray_gate_comp.png');