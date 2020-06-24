%WG_Plot_Test.m
clc;
clear all;
close all;

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
load('XY_quantum_gates.mat','-mat')

%Waveguides
n=5;
N=n-1;  %Photons? classically

%Trotter Iteration number
steps=1;

%Iteration time
t=pi/2;

%Create initial pure state
m=n;
theta=pi;    %relative phase
phi0=NOON_Wave(1,n,m,theta);

%% Exact decomposition
[Hij]=XY_Exchange_Hamiltonian(n);
%FockPrint(Hij)
s=50;
phi_i=zeros(length(phi0),s);
phi_i(:,1)=phi0;
ts=linspace(0,t,s);
for i=2:s
    phi_i(:,i)=expm(-1i*Hij*ts(i))*phi0;
end

indexes=2.^(0:n-1)+1;
phi_i_red=phi_i(indexes,:);

mode={'fock_bin',1,n};
ket_string=ket_fock_str(indexes,mode);
%ket_string=ket_fock_tex(indexes,mode); %-> for use with latex

%Plot waveguide development of ideal waveguide
WG_Plot2PColor( ket_string,phi_i_red,ts);