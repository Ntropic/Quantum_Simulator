%test_bintree.m
clc;
clear all;
close all;
fclose('all');

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened));

%Create binary tree
depth=6;
n=33;
file_name='tree';
mode=[1 0];

tree=bintree(depth,n);

if mode(1)==1
    path=curr_path(tree)
    path_rec=pointer2path(tree,6)

    %Plot the tree
    plot_obj=plot_tree(tree);
    fileID=fopen([file_name,'.tex'],'w');
    fprintf(fileID,plot_obj);
    fclose(fileID);

    tex2pdf2preview(file_name,1,400);
end
if mode(2)==1
    %Modify Layer
    [tree2 lay_state]=layer_up(tree);
    [tree2 lay_state]=layer_up(tree2);
    [tree2 lay_state]=layer_up(tree2);

    %Plot the modified tree
    plot_obj=plot_tree(tree2);
    fileID=fopen([file_name,'2.tex'],'w');
    fprintf(fileID,plot_obj);
    fclose(fileID);

    tex2pdf2preview([file_name,'2'],1,400);
end