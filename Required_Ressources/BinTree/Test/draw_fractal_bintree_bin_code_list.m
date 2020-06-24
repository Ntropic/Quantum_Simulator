%test_bintree_bin_code_list.m
clc;
clear all;
close all;
fclose('all');

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened));

%Create binary tree
depth=3;
n=1;
file_name='tree';
mode=[1 0];

tree=bintree(depth,n);

perm=all_num_algorithm(tree);
dec2bin(perm-1)