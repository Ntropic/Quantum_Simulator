%test_bintree_second_animation.m
clc;
clear all;
close all;
fclose('all');

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened));

mode=2;

if mode==1
    depth=6;
    for n=1:2^8-1

        tree=bintree(depth,n);

        perm=all_num_algorithm(tree);

        if n>1
            h=permutation_draw(perm,'lines',h);
        else
            h=permutation_draw(perm,'lines');
        end
        pause(0.3)
    end
elseif mode==2
    %Create binary tree
    for i=2:14
        depth=i;
        n=1;

        tree=bintree(depth,n);

        perm=all_num_algorithm(tree);

        if i>2
            h=permutation_draw(perm,'lines',h);
        else
            h=permutation_draw(perm,'lines');
        end
        a=get(h,'children');
        title(a(1),['d=' num2str(permutation_dimension(perm))]);
        pause(0.3)
    end
end