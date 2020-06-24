%test_markovian_noise.m
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

printing=0;

%% 1st: Single Qubit Decoherence
%Model details: rho0=[1-rho11,rho01;rho01*,rho11]
%               ---> [1-rho11*e^(-t/T1),rho01*e^(-t/(2T1))*e^(-(t/T2)^(1+alpha)))

%n=1;
n=4;

%Error model
t_step=0.01;
T1=10;
T2=1;

%[ K,p,ind,chan ] = Pauli_Twirling( n,t_step,T1,T2 );
[ K,chan ] = Kraus_Twirling( n,t_step,T1,T2 );

[ K2 ] = Kraus_Twirling( n/2,t_step,T1,T2 );
for i=1:length(K2)
    K3{i}=kron(K2{i},diag(ones(2^(n/2),1)));
    K4{i}=kron(diag(ones(2^(n/2),1)),K2{i});
end
%% Decoherence
%wave=[0,1];%1/sqrt(2)*[1,-1i];
wave=fock2vec(2,2,[2 0;0 2],[1 1i]/sqrt(2));

rho=sparse(Wave2Density(wave,1));

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

rho3=rho;
t1=tic();
for steps=1:100
    rho2=rho;
    rho=rho*0;
    %rho=rho*p(1);
    for i=1:length(K)
        rho=rho+K{i}*rho2*K{i}';
        %rho=rho+p(i+1)*K{i}*rho2*K{i}';
    end
    
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
toc(t1)
%[~,~,M]=Sparse2Square(rho,{'fock_full',2,4});
%fprintf(M)

t2=tic();
for steps=1:100
    rho4=rho3;
    rho3=rho3*0;
    %rho=rho*p(1);
    for i=1:length(K3)
        rho3=rho3+K3{i}*rho4*K3{i}';
    end
    
    rho4=rho3;
    rho3=rho3*0;
    
    for i=1:length(K4)
        rho3=rho3+K4{i}*rho4*K4{i}';
    end
end
toc(t2)

x=rho3-rho;
max(abs(x(:)))