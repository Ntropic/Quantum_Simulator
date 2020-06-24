%WG_Plot_Test_Compare_Depth.m
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
n=8;
N=n-1;  %Photons? classically

depth=7;
s_e=1;
d_e=1;

%Trotter Iteration number
steps=3;

%Iteration time
t=pi/4;

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

%% Approximate decomposition
%% Prepare gates
Ni=N:-1:1;
in=1:N;
Jx=sqrt(in).*sqrt(Ni);

%Gate A
sizle=floor(n/2)*2;
anc_size=0;
nameA=['A_n_' num2str(n)];
circ_string=['A_{XY}^{n=' num2str(n) '}'];
gates_indexes={};
param_steps={};
gates_names={};
for k=1:2:2*floor(n/2)-1
    gates_indexes={gates_indexes{:},[k,k+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
%A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

%Gate B
sizle=floor((n-1)/2)*2;
anc_size=0;
nameB=['B_n_' num2str(n)];
circ_string=['B_{XY}^{n=' num2str(n) '}'];
gates_indexes={};
param_steps={};
gates_names={};
for k=1:2:2*floor((n-1)/2)-1
    gates_indexes={gates_indexes{:},[k,k+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k+1)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
%B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);
        
%Create Trotter steps
name=['XY_n_' num2str(n)];
circ_string=['XY_{\text{approx}}^{n=' num2str(n) '}'];
circ_string_step=['XY_{\text{step}}^{n=' num2str(n) '}'];
gates_names={nameA,nameB};
gates_indexes={[1:A.size],1+[1:B.size]};
other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
%This command creates a Trotter gate for all the steps -> Great for
[gate,step_gate trotter_name trotter_step_name]=Trotter(t,steps,name,circ_string,circ_string_step,n,0,gates_names,gates_indexes);

%Do Trotter steps
other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

ind=Gate_Index_by_Name(trotter_name,elem_gates,other_gates);
[phi_a,ta]=DepthExpansion2Matrix_Return_Wave(phi0,other_gates(ind(2)),depth,s_e,d_e);
phi_a_red=phi_a(indexes,:);


%% Plotting 
mode={'fock_bin',1,n};
ket_string=ket_fock_str(indexes,mode);
%indexes2
mode2={'fock',n-1,2};
l=length(dec2bin(length(indexes)-1));
for i=1:length(indexes)
    d=[i-1,length(indexes)-i+1];
    indexes2(i)=bin2dec([dec2bin(d(1),l) dec2bin(d(2),l)]);
end
ket_string2=ket_fock_str(indexes2,mode2);
%ket_string=ket_fock_tex(indexes,mode); %-> for use with latex

%Plot waveguide development of ideal waveguide
WG_Plot2PColor_Compare( ket_string,ket_string2,phi_i_red,phi_a_red,ts,ta,0,1,50);