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
n_m=2:6;
n_h=2:8;
%Trotter Iteration number
steps=100;
s_h=10;

%Iteration time
t=pi/4;
h=figure();
ax=axes(h);

%% ------------------------------------------------------------------------------------------------------------------------
sym_name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];
trot_name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_XY'];
higher_name=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(30) '_r_1_t_' num2str(t) '_XY'];
higher_name_1=['higher_compare_n_' num2str(min(n_h)) '_to_' num2str(max(n_h)) '_s_' num2str(s_h) '_r_1_r_l_2_t_' num2str(t) '_XY'];
higher_name_r0=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(30) '_r_0_t_' num2str(t) '_XY'];
higher_name_r0_1=['higher_compare_n_' num2str(min(n_h)) '_to_' num2str(max(n_h)) '_s_' num2str(s_h) '_r_0_r_l_2_t_' num2str(t) '_XY'];
higher_name_r10_1=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(s_h) '_r_1_0_t_' num2str(t) '_XY'];
higher_name_r01_1=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(s_h) '_r_0_1_t_' num2str(t) '_XY'];
higher_name_r2=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(15) '_r_2_t_' num2str(t) '_XY'];
higher_name_r22=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(s_h) '_r_2_2_t_' num2str(t) '_XY'];
higher_name_r3=['higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(15) '_r_3_t_' num2str(t) '_XY'];

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

if exist(['../Higher_Order_Trotter/' higher_name_1 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_1 '.mat']);
    F_av31_ind=a.F_av_ind;
    gate_num31=a.gate_num;
    s_q31=a.s_q;
    fprintf('Found Higher Order Trotter calculations (l_r=2) \n')
else
    error('Higher not found (l_r=2)')
end

if exist(['../Higher_Order_Trotter/' higher_name_r0 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r0 '.mat']);
    F_av4_ind=a.F_av_ind;
    gate_num4=a.gate_num;
    s_q4=a.s_q;
    fprintf('Found Higher Order Trotter (r=0) calculations\n')
else
    error('Higher (r=0) not found')
end

if exist(['../Higher_Order_Trotter/' higher_name_r0_1 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r0_1 '.mat']);
    F_av5_ind=a.F_av_ind;
    gate_num5=a.gate_num;
    s_q5=a.s_q;
    fprintf('Found Higher Order Trotter (r=0 r_l=2) calculations\n')
else
    error('Higher (r=0 r_l=2) not found')
end

if exist(['../Higher_Order_Trotter/' higher_name_r10_1 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r10_1 '.mat']);
    F_av6_ind=a.F_av_ind;
    gate_num6=a.gate_num;
    s_q6=a.s_q;
    fprintf('Found Higher Order Trotter (r=1 0 r_l=2) calculations\n')
else
    error('Higher (r=1 0 r_l=2) not found')
end

if exist(['../Higher_Order_Trotter/' higher_name_r01_1 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r01_1 '.mat']);
    F_av61_ind=a.F_av_ind;
    gate_num61=a.gate_num;
    s_q61=a.s_q;
    fprintf('Found Higher Order Trotter (r=0 1 r_l=2) calculations\n')
else
    error('Higher (r=0 1 r_l=2) not found')
end

if exist(['../Higher_Order_Trotter/' higher_name_r2 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r2 '.mat']);
    F_av7_ind=a.F_av_ind;
    gate_num7=a.gate_num;
    s_q7=a.s_q;
    fprintf('Found Higher Order Trotter (r=2 ) calculations\n')
else
    error('Higher (r=2 ) not found')
end


if exist(['../Higher_Order_Trotter/' higher_name_r22 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r22 '.mat']);
    F_av71_ind=a.F_av_ind;
    gate_num71=a.gate_num;
    s_q71=a.s_q;
    fprintf('Found Higher Order Trotter (r=2 2) calculations\n')
else
    error('Higher (r=2 2) not found')
end

if exist(['../Higher_Order_Trotter/' higher_name_r3 '.mat'])
    a=load(['../Higher_Order_Trotter/' higher_name_r3 '.mat']);
    F_av8_ind=a.F_av_ind;
    gate_num8=a.gate_num;
    s_q8=a.s_q;
    fprintf('Found Higher Order Trotter (r=3) calculations\n')
else
    error('Higher (r=3) not found')
end

%% Analyze -------------------------------------------------------------------------------------------------------------
n=6;
%Plot 1-Fidelity over Steps and Photons
f_1=1-F_av_ind(n-1,:)';
f_2=1-F_av2_ind(n-1,:)';
f_3=1-F_av3_ind(n-2,:)';
f_31=1-F_av31_ind(n-2,:)';
f_4=1-F_av4_ind(n-2,:)';
f_5=1-F_av5_ind(n-2,:)';
f_6=1-F_av6_ind(n-2,:)';
f_61=1-F_av61_ind(n-2,:)';
f_7=1-F_av7_ind(n-2,:)';
f_71=1-F_av71_ind(n-2,:)';
f_8=1-F_av8_ind(n-2,:)';
[x1,n1]=Order_Of_Convergence(f_1,1:steps);
[x2,n2]=Order_Of_Convergence(f_2,1:steps);
[x3,n3]=Order_Of_Convergence(f_3,1:30);
[x31,n31]=Order_Of_Convergence(f_31,1:s_h);
[x4,n4]=Order_Of_Convergence(f_4,1:30);
[x5,n5]=Order_Of_Convergence(f_5,1:s_h);
[x6,n6]=Order_Of_Convergence(f_6,1:s_h);
[x61,n61]=Order_Of_Convergence(f_61,1:s_h);
[x7,n7]=Order_Of_Convergence(f_7,1:s_h);
[x71,n71]=Order_Of_Convergence(f_71,1:s_h);
[x8,n8]=Order_Of_Convergence(f_8,1:s_h);
plot(n1,x1,'-.r')
hold on
plot(n2,x2,'r')

plot(n4,x4,'--k')
plot(n3,x3,'--b')
plot(n7,x7,'--m')
plot(n8,x8,'--c')

plot(n5,x5,'k')
plot(n31,x31,'b')
plot(n71,x71,'m')

plot(n6,x6,'-.g')
plot(n61,x61,'g')


xlabel('\# Trotter steps')
ylabel('Order of Convergence')
yl=get(gca,'ylim');
axis([2 10 0 15]);
%legend({'$\hat{T}$','$\hat{S}_{1}$','$\hat{S}^{(3)}_{3}$','$\hat{S}^{(5)}_{3}$','$\hat{S}^{(7)}_{3}$','$\hat{S}^{(9)}_{3}$','$\hat{S}^{(3,3)}_{5}$','$\hat{S}^{(5,5)}_{5}$','$\hat{S}^{(7,7)}_{5}$','$\hat{S}^{(3,5)}_{5}$','$\hat{S}^{(5,3)}_{5}$'},'Interpreter','latex')
drawnow;

matlab2tikz(['conv_xy_all_comp_n_' num2str(n) '_t_' num2str(t) '.tex'],'parseStrings',false,'figurehandle',h,'standalone',true)

%% Plot the error over gate numbers
%Plot 1-Fidelity over Steps and Photons
gn_1=gate_num(n-1,:)';
gn_2=gate_num2(n-1,:)';
gn_3=gate_num3(n-2,:)';
gn_31=gate_num31(n-2,:)';
gn_4=gate_num4(n-2,:)';
gn_5=gate_num5(n-2,:)';
gn_6=gate_num6(n-2,:)';
gn_61=gate_num61(n-2,:)';
gn_7=gate_num7(n-2,:)';
gn_71=gate_num71(n-2,:)';
gn_8=gate_num8(n-2,:)';

h2=figure();
hold off;
plot(gn_1,f_1,'-.r')
hold on
plot(gn_2,f_2,'r')

plot(gn_4,f_4,'k--')
plot(gn_3,f_3,'b--')
plot(gn_7,f_7,'--m')
plot(gn_8,f_8,'--c')

plot(gn_5,f_5,'k')
plot(gn_31,f_31,'b')
plot(gn_71,f_71,'m')

plot(gn_6,f_6,'-.g')
plot(gn_61,f_61,'g')
legend({'$\hat{T}$','$\hat{S}_{1}$','$\hat{S}^{(3)}_{3}$','$\hat{S}^{(5)}_{3}$','$\hat{S}^{(7)}_{3}$','$\hat{S}^{(9)}_{3}$','$\hat{S}^{(3,3)}_{5}$','$\hat{S}^{(5,5)}_{5}$','$\hat{S}^{(7,7)}_{5}$','$\hat{S}^{(3,5)}_{5}$','$\hat{S}^{(5,3)}_{5}$'},'Interpreter','latex')

xlabel('\# Gates')
set(gca, 'YScale', 'log')
ylabel('1-F')
yl=get(gca,'ylim');
set(gca,'xlim',[0 9000]);
drawnow;

drawnow;
matlab2tikz(['error_gates_xy_all_comp_n_' num2str(n) '_t_' num2str(t) '.tex'],'parseStrings',false,'figurehandle',h2,'standalone',true)


% matlab2tikz('xy_gate_comp.tex','figurehandle',h2)
% pause() 
% A=getframe(h);
% B=getframe(h2);
% imwrite(A.cdata,'higher_xy_trotter_comp.png');
% imwrite(B.cdata,'higher_xy_gate_comp.png');