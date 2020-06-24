%test_bintree_animation.m
clc;
clear all;
close all;
fclose('all');

mode=1;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened));

%Create binary tree
if mode==1
    for depth=2
        for n=1:floor(2^depth)
            file_name='tree';
            gif_name=['algorithm_' num2str(depth) '_' num2str(n) '.gif'];

            tree=bintree(depth,n);

            perm=all_num_algorithm_gif(tree,gif_name)
        end
    end
else
    for depth=2
        %for n=1:floor(2^depth)

            tree=bintree(depth,1);

            perm=all_num_algorithm_tex(tree);
        %end
    end
end