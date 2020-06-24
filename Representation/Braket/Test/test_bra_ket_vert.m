%test_bra_ket_vert.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))


ket_vert=vertfock_str(1:128,{'fock unit dense',2,3,1});
fprintf(ket_vert{128})
fprintf('\n')

ket_vert=vertfock_str(1:128,{'fock_bin unit dense',2,3,1});
fprintf(ket_vert{128})
fprintf('\n')


ket=ket_fock_str(1:128,{'fock unit dense',2,3,1});
fprintf(ket{128})
fprintf('\n')

ket=ket_fock_str(1:128,{'fock_bin unit dense',2,3,1});
fprintf(ket{128})
fprintf('\n')

bra=bra_fock_str(1:128,{'fock unit',2,3,1});
fprintf(bra{128})
fprintf('\n')

bra=bra_fock_str(1:128,{'fock_bin unit',2,3,1});
fprintf(bra{128})
fprintf('\n')