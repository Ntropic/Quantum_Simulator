%test_extrapolation.m
clc;
clear all;
close all;

%Add paths
p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-4});
addpath(genpath(shortened));

%Load simulate data
load('mwg_trot_compare_fid_n_1_to_7_s_10_t_0.7854_N_4_pauli_ladder.mat');
load('mwg_trot_compare_n_1_to_7_s_10_t_0.7854_N_4_pauli_ladder.mat');

%Extrapolate convergence order from simulation data
wg=3;
ph=4;

f=1-F_av_ind(wg,ph,:);
o=Order_Of_Convergence(f);

[fitresult,gof]=Convergence_Extrapolation_Regression(o,1:length(o),0);
fita=fitresult.a;
fitb=fitresult.b;
fitc=fitresult.c;

%%Check if this reproduces the correct error rates
plot(f(:),'r');
f_fun=@(m,n,fn,c) fn*(n./m).^c;
f2=f_fun(7:20,6,f(6),fitc);
hold on;
plot(7:20,f2,'b')