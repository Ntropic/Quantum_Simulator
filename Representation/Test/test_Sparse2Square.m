%test_Sparse2Square.m
clc;
clear all;
close all;

testmode=3;


p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-2});
addpath(genpath(shortened));

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

%Test CreateControlledCompleteSet_Cn_U2 for n
%Create/Load Controlled Unitaries
n=3;
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);
S=PU_gates(11).matrix;

if testmode==1
    [~,~,M_mode]=Sparse2Square( S,'none');
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,{'fock',1,n});
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,'fock_dense');
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square( S,'index');
    fprintf(M_mode);
elseif testmode==2
    global alpha beta delta theta subs_comp
    alpha=sym('alpha');
    beta=sym('beta');
    delta=sym('delta');
    theta=sym('theta');
    assume(alpha,'real');
    assume(beta,'real');
    assume(delta,'real');
    assume(theta,'real');
    subs_comp={(exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,
               (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2};
    subs_list={'U11','U12',
               'U21','U22'};
    
    [~,~,M_mode]=Sparse2Square( S,'none_def_subs');
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,{'fock',1,n},'substitution_list',subs_list,'substitution_compare',subs_comp);
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,{'fock_dense',1,n},'substitution_list',subs_list,'substitution_compare',subs_comp);
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square( S,'index','substitution_list',subs_list,'substitution_compare',subs_comp);
    fprintf(M_mode);
elseif testmode==3
    [~,~,M_mode]=Sparse2Square( S,'none_def_subs_var_width');
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,{'fock_def_subs_var_width',2,2});
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square(S,'fock_dense_def_subs_var_width');
    fprintf(M_mode);
    [~,~,M_mode]=Sparse2Square( S,'index_def_subs_var_width');
    fprintf(M_mode);
end