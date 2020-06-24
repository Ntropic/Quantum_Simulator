%test_2_markovian_noise.m
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

K=Efficient_Kraus_Twirling( n,t_step,T1,T2 );

%% Decoherence
%wave=[0,1];%1/sqrt(2)*[1,-1i];
wave=fock2vec(1,3,[1 0  0;0 0 1],[1 1i]/sqrt(2));

rho=sparse(Wave2Density(wave,1));

if printing==1
    b=bar3(abs(rho));

    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end

    axis tight;
    x=xlim;
    y=ylim;
    z=[0 1];%zlim;
    xlim(x);
    ylim(y);
    zlim(z);
    asp=daspect;
    asp(2)=asp(1);
    daspect(asp);
    drawnow()
    %pause()
    %matlab2tikz('start_q.tex')
end

rho3=rho;
for steps=1:100
    rho=Kraus_Twirl_Step( rho,n,K );
    
    if printing==1
        b=bar3(abs(rho));

        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end

        xlim(x);
        ylim(y);
        zlim(z);
        daspect(asp);
        drawnow()
    end
end

%matlab2tikz('end_q.tex')
%[~,~,M]=Sparse2Square(rho,{'fock_full',2,4});
%fprintf(M)