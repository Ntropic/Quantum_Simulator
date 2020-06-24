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
%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
load('../clifford_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

%Photon number
n_m=1:4;

%Trotter Iteration number
steps=5;
r_h=0;
r_h2=1;

%Iteration time
t=pi/4;

h=figure('Position', [10 10 1000 700]);
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
name=['comp_sym_trot_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_ladder'];
sym_name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_ladder'];
trot_name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_ladder'];
higher_name=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_r_' num2str(r_h) '_t_' num2str(t) '_pauli_ladder'];
higher_name2=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_r_' num2str(r_h2) '_t_' num2str(t) '_pauli_ladder'];
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
n=4;
%Plot 1-Fidelity over Steps and Photons
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
f_3=1-F_av3_ind(2:end,:)';
f_4=1-F_av4_ind(2:end,:)';
[x1,n1]=Order_Of_Convergence(f_1,1:steps);
[x2,n2]=Order_Of_Convergence(f_2,1:steps);
[x3,n3]=Order_Of_Convergence(f_3,1:steps);
[x4,n4]=Order_Of_Convergence(f_4,1:steps);
plot(n1,x1,'Color',[1 0.5 0])
hold on
plot(n2,x2,'g')
plot(n3,x3,'b')
plot(n4,x4,'r')
legend({'Simple','Symmetric','Tripling','Quintupling'})
xlabel('# Trotter steps')
ylabel('Order')
axis([2 5 0 9]);
drawnow;

drawnow;
matlab2tikz(['conv_ladder_all_comp_n_' num2str(n) '.tex'],'figurehandle',h)
% matlab2tikz('xy_gate_comp.tex','figurehandle',h2)
% pause() 
% A=getframe(h);
% B=getframe(h2);
% imwrite(A.cdata,'higher_xy_trotter_comp.png');
% imwrite(B.cdata,'higher_xy_gate_comp.png');