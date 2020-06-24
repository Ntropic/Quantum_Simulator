%test_speed_symbolic_Sparse2Square.m
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

%Test CreateControlledCompleteSet_Cn_U2 for n
n=4;

alpha=sym('alpha');
beta=sym('beta');
delta=sym('delta');
theta=sym('theta');
assume(alpha,'real');
assume(beta,'real');
assume(delta,'real');
assume(theta,'real');
subs_comp={(exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2;
           (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2};


%Generate Gates
t=tic();
[PU_gates,Cn]=CreateCompleteSet_Cn_U2(n,elem_gates,comp_gates,[1,1,0]);
t2=toc(t);
fprintf('Created %1d gates in %1.2f s.\n',size(PU_gates,1)*size(PU_gates,2),t2)

%Print Gates
for i=1%:size(PU_gates,1)
    for j=1%:size(PU_gates,2)
        %Two methods
        M=PU_gates(i,j).matrix;
        M_str=cell(size(M));
        
        %First
        t=tic();
        M_prior=char(M);
        M_prior=M_prior(10:end-2);
        M_vec=strsplit(M_prior,'], [');
        M_str=regexp(M_vec,', ','split');
        M_str=vertcat(M_str{:});
        [ind1 ind2]=find(strcmp(M_str,'1')+strcmp(M_str,'0')==0);
        time1=toc(t)
        
        %Second
        t=tic();
        for k=1:size(M,1)
            for l=1:size(M,2)
                M_str{k,l}=char(M(k,l));
            end
        end
        time2=toc(t)
        
    end
end
speedup=time2/time1
