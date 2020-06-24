%test_Sparse2Tex.m
clc;
clear all;
close all;

testmode=1;


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
    [filename]=Sparse2Tex(S,'Cn_test','none');
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test',{'fock',1,n});
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test','fock_dense');
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test','index');
    tex2pdf2preview(filename,1,600);
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
    subs_comp=[(exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,
               (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2];
    subs_list={'\hat{U}_{11}','\hat{U}_{12}',
               '\hat{U}_{21}','\hat{U}_{22}'};
    
    [filename]=Sparse2Tex(S,'Cn_test1','none_def_subs');
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test2',{'fock',1,n},'substitution_list',subs_list,'substitution_compare',subs_comp);
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test3',{'fock_dense',1,n},'substitution_list',subs_list,'substitution_compare',subs_comp);
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test4','index','substitution_list',subs_list,'substitution_compare',subs_comp);
    tex2pdf2preview(filename,1,600);
elseif testmode==3
    [filename]=Sparse2Tex(S,'Cn_test1','none_def_subs_var_width');
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test2',{'fock_def_subs_var_width',1,n});
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test3','fock_dense_def_subs_var_width');
    tex2pdf2preview(filename,1,600);
    [filename]=Sparse2Tex(S,'Cn_test4','index_def_subs_var_width');
   	tex2pdf2preview(filename,1,600);
end