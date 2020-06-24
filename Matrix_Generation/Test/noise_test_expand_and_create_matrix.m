%noise_test_expand_and_create_matrix.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

plotting=0;

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')


%% Create List with subgates & indexes
Toffoli11=Gate_by_Name('Toffoli11',elem_gates,comp_gates);

%Create a gate expansion into the matrix multiplication steps -> Future:
%Add Embedded Gates matrices
[gates_matrices gates_expanded gates_indexes_expanded N N_anc]=Gates_Expand(elem_gates,comp_gates,Toffoli11);

%% Noise Model

n=N_anc;

%Error model
t_step=0.0001;
T1=5;
T2=1;

phi=0;
E=0.0003;    %Error rate
delta=sqrt(10/3*E); 
E1=5/4*E;
%E=1-F -> E=2/5*E1+3/20*delta^2 -> with F=0.997 -> E=0.003 

K=Efficient_Kraus_Twirling( n,t_step,T1,T2 );
[K2,p2]=CZ_Twirling(phi,delta,E1);

%% Density Function
wave=fock2vec(1,3,[1 0 0;0 1 0;0 0 1],[1 1 1]/sqrt(3));
rho=full(sparse(Wave2Density(wave,1)));

rho2=double(Toffoli11.matrix*rho*Toffoli11.matrix');

for steps=1:length(gates_expanded)
    U=Gate_by_Name(gates_expanded{steps},elem_gates,comp_gates);
    twirl=U.twirling;
    index=gates_indexes_expanded{steps};
    U_emb=Embed_Gate(double(gates_matrices{steps}),index,n);
    
    rho=U_emb*rho*U_emb';
    if any(strfind(twirl,'cz'))
        rho=CZ_Twirl_Step( rho, K2, p2, index );
    end
    if any(strfind(twirl,'0'))~=1
        rho=Kraus_Twirl_Step( rho,n,K );
    end
end

tr_rho=trace(rho)
tr_rho2=trace(rho2)

%Fidelity between 2 Density Matrices
mixmax=max(max(abs(rho2-rho)))
F=(trace(sqrtm(sqrtm(rho2)*rho*sqrtm(rho2))))^2 %Numerically instable -> use other one instead
F=(trace(sqrtm(sqrtm(rho)*rho2*sqrtm(rho))))^2