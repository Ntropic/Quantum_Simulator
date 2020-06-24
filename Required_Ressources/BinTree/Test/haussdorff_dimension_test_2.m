%haussdorff_dimension_test_2.m
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


iter=50;
perm_len=1;
len=1;
dim(1)=1;
for i=2:iter
    len=2^i-1;
    perm_len=perm_len*2+2^(i-1);
    dim(i)=log(perm_len)/log(len);
end

perm_len=1;
len=1;
dim2(1)=1;
for i=2:iter
    len=2^i-1;
    perm_len=perm_len*2+2^(i-1)+len;
    dim2(i)=log(perm_len)/log(len);
end

perm_len=sqrt(2);
len=1;
dim3(1)=1;
for i=2:iter
    len=2^i-1;
    perm_len=perm_len*2+sqrt(2^((i-1)*2)+1);
    dim3(i)=log(perm_len)/log(len);
end

plot(1:iter,dim,'r')
hold on;
plot(1:iter,dim2,'b')
plot(1:iter,dim3,'g')