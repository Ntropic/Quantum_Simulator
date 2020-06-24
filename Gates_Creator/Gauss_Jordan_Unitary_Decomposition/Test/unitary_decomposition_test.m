%unitary_angle_matrix_test.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened))

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

alpha=sym('alpha','real');
beta=sym('beta','real');
theta=sym('theta','real');
delta=sym('delta','real');

U_gates.size=1;
U_gates.aux_size=0;
U_gates.step_num=4;
U_gates.steps.index={1,1,1,1};
U_gates.steps.gates={'EZ','EY','EZ','P01'};
U_gates.steps.param={[{sym('phi'),sym('beta','real')}],[{sym('phi'),sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),sym('delta','real')}]};

U_gates.matrix=Gate2Matrix(elem_gates,[],U_gates,2);
U_gates.matrix


U_gates2.size=1;
U_gates2.aux_size=0;
U_gates2.step_num=8;
U_gates2.steps.index={1,1,1,1,1,1,1,1};
U_gates2.steps.gates={'EZ','X','EZ','EY','X','EY','EZ','P01'};
U_gates2.steps.param={[{sym('phi'),1/2*(sym('beta','real')-sym('alpha','real'))}],[],[{sym('phi'),-(sym('beta','real')+sym('alpha','real'))/2}],[{sym('phi'),-1/2*sym('theta','real')}],[],[{sym('phi'),1/2*sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),sym('delta','real')}]};

U_gates2.matrix=Gate2Matrix(elem_gates,[],U_gates2,2);
U_gates2.matrix


U_gates3.size=2;
U_gates3.aux_size=0;
U_gates3.step_num=8;
U_gates3.steps.index={1,[1,2],1,1,[1,2],1,1,2};
U_gates3.steps.gates={'EZ','CNOT','EZ','EY','CNOT','EY','EZ','P1'};
U_gates3.steps.param={[{sym('phi'),1/2*(sym('beta','real')-sym('alpha','real'))}],[],[{sym('phi'),-(sym('beta','real')+sym('alpha','real'))/2}],[{sym('phi'),-1/2*sym('theta','real')}],[],[{sym('phi'),1/2*sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),sym('delta','real')}]};

U_gates3.matrix=Gate2Matrix(elem_gates,[],U_gates3,2);
U_gates3.matrix

a=2*pi*rand(1);
b=2*pi*rand(1);
d=pi*rand(1);
t=2*pi*rand(1);

U1=subs(U_gates.matrix,alpha,a);
U1=subs(U1,beta,b);
U1=subs(U1,theta,t);
U1=subs(U1,delta,d);
U1_test=double(U1)

U2=subs(U_gates2.matrix,alpha,a);
U2=subs(U2,beta,b);
U2=subs(U2,theta,t);
U2=subs(U2,delta,d);
U2_test=double(U2)


U3=subs(U_gates3.matrix,alpha,a);
U3=subs(U3,beta,b);
U3=subs(U3,theta,t);
U3=subs(U3,delta,d);
U3_test=double(U3)