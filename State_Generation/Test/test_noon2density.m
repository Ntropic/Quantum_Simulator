%test_noon2density.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

wave=fock2vec(2,4,[2 0 0 0; 0 0 0 2 ],[1 1i]/sqrt(2))

rho=sparse(Wave2Density(wave,1))

[~,~,M]=Sparse2Square(rho,{'fock',2,4});
fprintf(M)