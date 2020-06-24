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

for i=1:100
%Create unitary matrix
alpha=2*pi*rand(1);
beta=2*pi*rand(1);
% phi1=mod(-1/2*(alpha+beta)+2*pi,2*pi);
% phi2=mod(-1/2*(alpha-beta)+2*pi,2*pi);
delta=pi*rand(1);
theta=pi*rand(1);

U=[cos(theta/2)*exp(1i*(delta-alpha/2-beta/2)), -sin(theta/2)*exp(1i*(delta-alpha/2+beta/2)) ; sin(theta/2)*exp(1i*(delta+alpha/2-beta/2)) , cos(theta/2)*exp(1i*(delta+alpha/2+beta/2))];
%U=exp(1i*delta)*[cos(theta/2)*exp(1i*p1), -sin(theta/2)*exp(1i*p2) ; sin(theta/2)*exp(-1i*p2) , cos(theta/2)*exp(-1i*p1)]

[alpha2,beta2,delta2,theta2]=Unitary2Angles(U);
%[alpha2,beta2,delta2,theta2]

U2=[cos(theta2/2)*exp(1i*(delta2-alpha2/2-beta2/2)), -sin(theta2/2)*exp(1i*(delta2-alpha2/2+beta2/2)) ; sin(theta2/2)*exp(1i*(delta2+alpha2/2-beta2/2)) , cos(theta2/2)*exp(1i*(delta2+alpha2/2+beta2/2))];

U-U2

%U_gates=Angles2Unitary( elem_gates,comp_gates,alpha2,beta2,delta2,theta2 );

%U-U_gates.matrix
end