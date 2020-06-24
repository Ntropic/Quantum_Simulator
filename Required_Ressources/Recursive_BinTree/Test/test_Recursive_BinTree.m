%test_Recursive_BinTree.m
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
test=3;
m0=1;
N=5;

if test==0
    m=m0;
    perm=Rec_BinTree(N,m);

    %Compare
    tree=bintree(N,m); 
    perm_b=all_num_algorithm(tree);

    plot(m:2^N,perm,'r')
    hold on;
    plot(m:2^N,perm_b,'b')
    hold off;
    legend({'Recursive','Iterative'})
elseif test==1
    for m=m0:2^N
        perm=Rec_BinTree(N,m);

        %Compare
        tree=bintree(N,m); 
        perm_b=all_num_algorithm(tree);

        plot(m:2^N,perm,'r')
        hold on;
        plot(m:2^N,perm_b,'b')
        hold off;
        legend({'Recursive','Iterative'})
        pause(0.3);
    end
elseif test==2
    m=m0;
    t1=tic();
    perm=Rec_BinTree(N,m);
    t_rec=toc(t1)
    
    %Compare
    t2=tic();
    tree=bintree(N,m); 
    perm_b=all_num_algorithm(tree);
    t_tree=toc(t2)
    
    plot(m:2^N,perm,'r')
    hold on;
    plot(m:2^N,perm_b,'b')
    hold off;
    legend({'Recursive','Iterative'})
elseif test==3
    m=m0;
    N0=12;
    for N=1:N0
        t1=tic();
        perm=Rec_BinTree(N,m);
        t_rec(N)=toc(t1);

        %Compare
        t2=tic();
        tree=bintree(N,m); 
        perm_b=all_num_algorithm(tree);
        t_tree(N)=toc(t2);
    end
    
    semilogy(1:N0,t_rec,'r')
    hold on;
    semilogy(1:N0,t_tree,'b')
    hold off;
    legend({'Recursive','Iterative'})
end