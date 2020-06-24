%compare_test_Recursive_BinTree.m
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));


%% Test
m=1;
N=5;

    perm=Rec_BinTree(N,m);
    perm2=Rec_BinTree_2(N,m);

    %Compare
    tree=bintree(N,m); 
    perm_b=all_num_algorithm(tree);

    plot(m:2^N,perm,'r')
    hold on;
    plot(m:2^N,perm2,'b')
    plot(m:2^N,perm_b,'b')
    hold off;
    legend({'Recursive','New Recursive','Iterative'})
