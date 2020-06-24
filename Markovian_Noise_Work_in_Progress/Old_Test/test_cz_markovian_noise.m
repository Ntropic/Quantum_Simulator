%test_cz_markovian_noise.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

printing=1;

%% 1st: Single Qubit Decoherence
%Model details: rho0=[1-rho11,rho01;rho01*,rho11]
%               ---> [1-rho11*e^(-t/T1),rho01*e^(-t/(2T1))*e^(-(t/T2)^(1+alpha)))

%n=1;
n=3;

%Error model
t_step=0.01;
T1=1;
T2=1;

phi=0;
E=0.003;    %Error rate
delta=sqrt(10/3*E); 
E1=5/4*E;
%E=1-F -> E=2/5*E1+3/20*delta^2 -> with F=0.997 -> E=0.003 

K=Efficient_Kraus_Twirling( n,t_step,T1,T2 );
[K2,p2]=CZ_Twirling(phi,delta,E1);

CZ=Gate_by_Name( 'CZ',elem_gates,comp_gates );
CZ=CZ.matrix;
%% Decoherence
wave=fock2vec(1,3,[1 0 0;0 1 0],[1 1i]/sqrt(2));

rho=full(sparse(Wave2Density(wave,1)));
if printing==1
    fprintf([num2str(1),'/',num2str(100),'\n']);
    [F,bra]=DensityBar3( rho,'Coloring',true,'AbsValue',true);
    A=Bar3Plot2Tikz2Pdf2Png(F,bra,'Aux_Files','frame_title',1260);  %titleframe
    imwrite(A,'frame_title.png')
    A=Bar3Plot2Tikz2Pdf2Png(F,bra,'Aux_Files','frame',1260);
    v=VideoWriter(['Videos/CZ_Markovian_n_',num2str(n),'.avi'],'Uncompressed AVI');
    open(v)
    writeVideo(v,A);
elseif printing==2
    [F,bra]=DensityBar3( rho,'Coloring',true,'AbsValue',true);
end

for steps=1:100
    fprintf([num2str(steps+1),'/',num2str(100),'\n']);
    index=[mod(steps,2)+1,mod(steps+1,2)+1];
    rho=Embed_Gate(double(CZ),index,n)*rho*Embed_Gate(double(CZ),index,n)';
    rho=CZ_Twirl_Step( rho, K2, p2, index );
    rho=Kraus_Twirl_Step( rho,n,K );
    
    if printing==1
        F=DensityBar3( rho,'Coloring',true,'AbsValue',true,'Fig',F);
        A=Bar3Plot2Tikz2Pdf2Png(F,bra,'Aux_Files','frame',1260);
        writeVideo(v,A);
    elseif printing==2
        F=DensityBar3( rho,'Coloring',true,'AbsValue',true,'Fig',F);
    end
end
close(v);
%matlab2tikz('end_cz.tex','mathmode',true,'parseStrings',false)

%[~,~,M]=Sparse2Square(rho,{'fock_full',2,4});
%fprintf(M)