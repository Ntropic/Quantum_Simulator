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
n_h=2:8;
%Trotter Iteration number
steps=10;
s_h=5;
r_h=1;
r2_h=5;


dt=0.0001;
E=0.005; 

%Iteration time
t=pi/4;
h=figure();
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
name=['comp_sym_trot_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];
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
    error('Trot not found')
end

%% -------------------------------------------------------------------------------------------------------------------------

if exist(['../Simple_Trotter/' trot_name '.mat'])
    a=load(['../Simple_Trotter/' trot_name '.mat']);
    F_av_ind=a.F_av_ind;
    gate_num=a.gate_num;
    s_q=a.s_q;
    fprintf('Found Simple Trotter calculations\n')
else
    error('Sym not found')
end

%% -------------------------------------------------------------------------------------------------------------------------

if exist(['../Higher_Order_Trotter/' higher_name '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name '.mat']);
    F_av3_ind=a.F_av_ind;
    gate_num3=a.gate_num;
    s_q3=a.s_q;
    fprintf('Found Higher Order Trotter calculations\n')
else
    error('Higher not found')
end

%% Analyze -------------------------------------------------------------------------------------------------------------
n=6;
%Plot 1-Fidelity over Steps and Photons
f_1=1-F_av_ind(2:end,:)';
f_2=1-F_av2_ind(2:end,:)';
f_3=1-F_av3_ind(1:end,:)';
[x1,n1]=Order_Of_Convergence(f_1,1:steps);
[x2,n2]=Order_Of_Convergence(f_2,1:steps);
[x3,n3]=Order_Of_Convergence(f_3,1:s_h);
plot(n1,x1,'Color',[1 0.5 0])
hold on
plot(n2,x2,'g')
plot(2:10,[x3 ones(1,5)*8],'b')
legend({'Simple Trotter','Symmetric Trotter','Quintupling Trotter'})
xlabel('# Trotter steps')
ylabel('Order')
drawnow;

drawnow;
matlab2tikz(['conv_xy_all_comp_n_' num2str(n) '.tex'],'figurehandle',h)
% matlab2tikz('xy_gate_comp.tex','figurehandle',h2)
% pause() 
% A=getframe(h);
% B=getframe(h2);
% imwrite(A.cdata,'higher_xy_trotter_comp.png');
% imwrite(B.cdata,'higher_xy_gate_comp.png');